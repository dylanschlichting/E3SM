! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!

module ocn_submesoscale_eddies

   use mpas_pool_routines
   use mpas_derived_types
   use mpas_constants
   use mpas_threading
   use mpas_timer

   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables
   use ocn_mesh

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_submesoscale_compute_velocity, &
             ocn_submesoscale_init,             &
             ocn_submesoscale_finalize,         &
             ocn_MLE_add_to_transport_vel,      &
             ocn_submesoscale_update_filtered_fields

   real(kind=RKIND) :: tau, LfMin, Ce, dsMax, m_star, n_star, Cr
   real(kind=RKIND) :: BLD_growing_Tfilt, BLD_decaying_Tfilt
   real(kind=RKIND) :: MLD_growing_Tfilt, MLD_decaying_Tfilt
   real(kind=RKIND),dimension(:),allocatable :: time_scale
   
   ! Storage for time-filtered boundary layer and mixed layer depths
   ! These are updated each time step to implement direction-dependent running mean
   real(kind=RKIND),dimension(:),allocatable :: BLD_filtered           ! Running mean of BLD (fast)
   real(kind=RKIND),dimension(:),allocatable :: BLD_filtered_slow      ! Running mean of BLD (slow)
   real(kind=RKIND),dimension(:),allocatable :: turbVelScale2_filtered ! Running mean of tur. vel. scale
   logical :: filters_initialized = .false.
   logical :: bodner_nan_reported = .false.
   ! Conservative CFL-style safety factor for the Bodner-induxed transport clip.
   real(kind=RKIND), parameter :: bodner_transport_cfl = 0.25_RKIND
  !***********************************************************************

  contains

  !***********************************************************************
  !
  !  routine ocn_submesoscale_compute_velocity
  !
  !> \brief
  !>   Computes the overturning streamfunction and induced transport
  !>   velocity associated with mixed-layer submesoscale eddies.
  !>
  !> \details
  !>   This routine implements the Fox-Kemper et al. (2011) and
  !>   Bodner et al. (2023) mixed-layer eddy parameterizations.
  !>   For the Fox-Kemper et al. (2011) parameterization, the
  !>   submesoscale overturning streamfunction is
  !>
  !>      Psi_FK =
  !>         C_e * Delta_s / L_f *
  !>         H^2 * |grad b| /
  !>         sqrt(f^2 + tau^-2) * mu(z)
  !> 
  !>   where
  !>
  !>      C_e      : nondimensional efficiency coefficient
  !>      Delta_s  : resolved horizontal grid scale,
  !>                 taken as min(dcEdge, dsMax)
  !>      L_f      : subgrid frontal width scale
  !>      H        : mixed-layer depth
  !>      f        : Coriolis parameter
  !>      tau      : equatorial regularization timescale
  !>      grad b   : mixed-layer averaged horizontal buoyancy gradient
  !>      mu(z)    : nondimensional vertical structure function
  !>
  !>   The frontal width scale L_f is parameterized following
  !>   Fox-Kemper et al. (2011) using the mixed-layer stratification,
  !>   buoyancy gradient, and Coriolis frequency.
  !>
  !>   The shape function mu(z) smoothly tapers the overturning
  !>   circulation to zero beneath the mixed layer.
  !>
  !>   For the Bodner et al. (2023) parameterization, the overturning
  !>   streamfunction is
  !>
  !>      Psi_Bod =
  !>         [ C_r * Delta_s * |f| * h * H^2 * |grad b| /
  !>           (m_* u_*^3 + n_* w_*^3)^(2/3) ] * mu(z)
  !>
  !>   where
  !>
  !>      C_r      : rescaled efficiency coefficient,
  !>                 C_r = C_e / C_L
  !>      C_L      : arrest scale coefficient, assumed here
  !>                 to be 0.25 following Bodner et al. (2023)
  !>      h        : boundary layer depth from KPP
  !>      u_*      : friction velocity,
  !>                 u_* = sqrt(|tau_wind| / rho_0) 
  !>      w_*      : convective velocity scale,
  !>                 w_* = (-B0 * h)^(1/3)
  !>      B0       : surface buoyancy flux
  !>      m_*,n_*  : empirical weighting coefficients
  !>
  !>   In this formulation, the frontal width is no longer diagnosed
  !>   from mixed-layer stratification. Instead, the overturning strength
  !>   is regulated by a turbulence-based arrest scale that depends on
  !>   both wind-driven and convective turbulent forcing. We draw heavily 
  !>   from MOM6's implementation to enforce numerical stability. This includes
  !>   double filtering the boundary layer depth, turbulent velocity scale, 
  !>   and limiting the streamfunction magnitude as shown in their  
  !>   src/parameterizations/lateral/MOM_mixed_layer_restrat.F90
  !
  !***********************************************************************

  subroutine ocn_submesoscale_compute_velocity(forcingPool,diagnosticsPool,dt)

     !-----------------------------------------------------------------
     !
     ! local variables
     !
     !-----------------------------------------------------------------

     type(mpas_pool_type), intent(in) :: forcingPool, diagnosticsPool
     real(kind=RKIND), intent(in) :: dt
     integer :: k, nEdges, nCells, iCell, iEdge, cell1, cell2, indMLDedge

      real(kind=RKIND),dimension(:),allocatable :: streamFunction, zEdge
      real(kind=RKIND) :: bodner_edge_scale, bodner_vel, bodner_vel_max

      ! Shared / FK2011 variables
     real(kind=RKIND) :: zMLD, bvfML, hML, gradBuoyML, bvfAv, hAv, ds

      if (trim(config_submesoscale_scheme) == 'Bod2023') then

         if (.not. filters_initialized) then
            normalMLEvelocity(:,:) = 0.0_RKIND
            return
         end if

      end if

     nEdges = nEdgesHalo(2)

     allocate(streamFunction(nVertLevels+1))
     allocate(zEdge(nVertLevels+1))

     ! compute the mixed layer average buoyancy gradient and bvf
     ! This is thickness weighted to allow for non uniform grid spacing, 
     ! and centered on layer interfaces. The mixed layer goes from the 
     ! ocean top layer to mid-depth of layer(indMLDedge)
     ! NOTE, gradDensityEddy and BFV are defined at the top cell interface 
     ! and not valid for minLevelCell, so properties for the top layer 
     ! rely on the values at (minLevelCell+1).

     !$omp parallel
     !$omp do schedule(runtime) &
     !$omp private(cell1,cell2,bvfML,zEdge,streamFunction,gradBuoyML,hML,bvfAv, &
     !$omp         zMLD,ds,k,indMLDedge,hAv)
     do iEdge=1,nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)

        zEdge(:) = 0.0_RKIND
        streamFunction(:) = 0.0_RKIND
        bvfML = 0.0_RKIND
        gradBuoyML = 0.0_RKIND

        !the BFV is defined at cell interfaces in the vertical so there is an offset in index
        !this is the top half layer, assuming homogeneity of top layer
        bvfAv = sqrt(0.5_RKIND*(max(BruntVaisalaFreqTop(minLevelCell(cell1)+1,cell1),1.0E-20_RKIND) + &
                       max(BruntVaisalaFreqTop(minLevelCell(cell2)+1,cell2),1.0E-20_RKIND)))

        hML = 0.5_RKIND*layerThickEdgeMean(minLevelEdgeBot(iEdge),iEdge)

        bvfML = bvfML + hML*bvfAv

        if (minLevelEdgeTop(iEdge) .ge. 1) then
           gradBuoyML = hML*gradDensityEddy(minLevelEdgeTop(iEdge)+1,iEdge)
           bvfML = hML*bvfAv
        else
           cycle
        end if

        indMLDedge = min(indMLD(cell1),indMLD(cell2))
        do k = minLevelEdgeBot(iEdge)+1,indMLDedge
           hAv = 0.5_RKIND*(layerThickEdgeMean(k,iEdge) + layerThickEdgeMean(k-1,iEdge))
           bvfAv = sqrt(0.5_RKIND*(max(BruntVaisalaFreqTop(k,cell1),1.0E-20_RKIND) + &
                              max(BruntVaisalaFreqTop(k,cell2),1.0E-20_RKIND)))
           bvfML = bvfML + hAv*bvfAv
           hML = hML + hAv
           gradBuoyML = gradBuoyML + hAv*gradDensityEddy(k,iEdge)
        end do

        bvfML = bvfML / (1.0E-20_RKIND + hML)
        gradBuoyML = gradBuoyML / (1.0E-20_RKIND + hML)
        gradBuoyML = gravity*gradBuoyML/rho_sw !convert from density unit to buoyancy

        !compute depths and shape function
        do k = minLevelEdgeTop(iEdge)+1,maxLevelEdgeTop(iEdge)+1
           zEdge(k) = zEdge(k-1) - layerThickEdgeMean(k-1,iEdge)
        end do

        zMLD = 0.5_RKIND*(dThreshMLD(cell1)+dThreshMLD(cell2))
        ds = min(dcEdge(iEdge),dsMax)

        if (config_submesoscale_scheme == 'FK2011') then
            call compute_streamfunction_fk2011( &
                iEdge,          &
                zMLD,           &
                ds,             &
                gradBuoyML,     &
                bvfML,          &
                zEdge,          &
                streamFunction)
        else if (config_submesoscale_scheme == 'Bod2023') then
            call compute_streamfunction_Bod2023( &
                iEdge,                &
                cell1,                &
                cell2,                &
                zMLD,                 &
                ds,                   &
                gradBuoyML,           &
                zEdge,                &
                streamFunction,       & 
                forcingPool,          &
                diagnosticsPool,      &
                dt,                   &
                layerThickEdgeMean,   &
                gradDensityEddy)
            ! ------------------------------------------------------------------
            ! MOM6-style transport limiter for the Bodner overturning profile.
            !
            ! The Bodner streamfunction Psi_k is converted to a layer-wise
            ! diagnostic transport tendency through
            !
            !   w_k = -(Psi_k - Psi_{k+1}) / h_k
            !
            ! where h_k = layerThickEdgeMean(k,iEdge).
            !
            ! This quantity represents the implied vertical transport associated
            ! with the vertical divergence of the MLE bolus circulation, rather
            ! than a prognostic vertical velocity.
            !
            ! The limiter checks the CFL-like constraint
            !
            !   |w_k| <= c * h_k / dt
            !
            ! with c = bodner_transport_cfl = 0.25. Larger c makes the limiter
            ! less restrictive.
            !
            ! Example caps for c = 0.25:
            !   - h_k = 10 m, dt = 300 s  -> |w_k| <= 8.3e-3 m/s  -> 2.5 m/step
            !   - h_k = 20 m, dt = 900 s  -> |w_k| <= 5.6e-3 m/s  -> 5.0 m/step
            !   - h_k = 50 m, dt = 1800 s -> |w_k| <= 6.9e-3 m/s  -> 12.5 m/step
            !
            ! These are per-timestep limits on implied tracer transport, so the
            ! accumulated displacement can still be O(100 m/day) or larger if the
            ! limiter is active every step.
            !
            ! If any layer violates this condition, we compute a single edge-wide
            ! scale factor
            !
            !   s = min_k [ (c * h_k / dt) / |w_k| ]
            !
            ! and then rescale the whole streamfunction:
            !
            !   Psi_k <- s * Psi_k
            !
            ! so the vertical shape of the overturning is preserved, but the
            ! amplitude is reduced until the implied transport constraint is
            ! satisfied everywhere in the column.
            !
            ! This mirrors the MOM6 psi_mag clip concept: clip the column-wide
            ! overturning strength, not each layer independently.
            ! ------------------------------------------------------------------
            bodner_edge_scale = 1.0_RKIND
            do k = minLevelEdgeTop(iEdge), maxLevelEdgeTop(iEdge)
               bodner_vel = -(streamFunction(k) - streamFunction(k+1)) / &
                            layerThickEdgeMean(k,iEdge)
               bodner_vel_max = bodner_transport_cfl * &
                                layerThickEdgeMean(k,iEdge) / max(dt, 1.0E-20_RKIND)
               if (abs(bodner_vel) > bodner_vel_max) then
                  bodner_edge_scale = min(bodner_edge_scale, &
                                          bodner_vel_max / max(abs(bodner_vel), 1.0E-20_RKIND))
               end if
            end do

            if (bodner_edge_scale < 1.0_RKIND) then
               streamFunction(:) = bodner_edge_scale * streamFunction(:)
            end if
        end if

        ! integrate in vertical to get the velocity
        do k = minLevelEdgeTop(iEdge),maxLevelEdgeTop(iEdge)
           normalMLEvelocity(k,iEdge) = -(streamFunction(k) - streamFunction(k+1)) / &
                                          layerThickEdgeMean(k,iEdge)
        end do

    end do!iEdge loop
    !$omp end do
    !$omp end parallel

    deallocate(streamFunction)
    deallocate(zEdge)

  end subroutine ocn_submesoscale_compute_velocity

  !***********************************************************************

  subroutine compute_streamfunction_fk2011( &
      iEdge, zMLD, ds, gradBuoyML, bvfML, zEdge, streamFunction)

      integer, intent(in) :: iEdge
      real(kind=RKIND), intent(in) :: zMLD, ds, gradBuoyML, bvfML
      real(kind=RKIND), dimension(:), intent(in) :: zEdge
      real(kind=RKIND), dimension(:), intent(inout) :: streamFunction

      integer :: k
      real(kind=RKIND) :: mu, Lf

      !   C_e*Delta_S/L_f*H^2 int(bGrad_H)/sqrt(f^2+tau^-2)*mu(z)
      Lf = max(Lfmin, &
               abs(gradBuoyML)*zMLD / time_scale(iEdge)**2.0_RKIND, &
               bvfML*zMLD / time_scale(iEdge))

      do k = minLevelEdgeTop(iEdge)+1,maxLevelEdgeTop(iEdge)

         mu = compute_mle_shape(zEdge(k), zMLD)

         streamFunction(k) = &
            Ce*ds/Lf*zMLD**2.0_RKIND* &
            gradBuoyML/time_scale(iEdge)*mu

      end do

  end subroutine compute_streamfunction_fk2011

  !***********************************************************************

  subroutine compute_streamfunction_Bod2023( &
      iEdge, cell1, cell2, zMLD_inst, ds, gradBuoyML_inst, zEdge, streamFunction, &
      forcingPool, diagnosticsPool, dt, layerThickEdgeMean_arr, gradDensityEddy_arr)

      integer, intent(in) :: iEdge, cell1, cell2
      real(kind=RKIND), intent(in) :: zMLD_inst, ds, gradBuoyML_inst, dt
      real(kind=RKIND), dimension(:), intent(in) :: zEdge
      real(kind=RKIND), dimension(:,:), intent(in) :: layerThickEdgeMean_arr
      real(kind=RKIND), dimension(:,:), intent(in) :: gradDensityEddy_arr
      real(kind=RKIND), dimension(:), intent(inout) :: streamFunction
      type(mpas_pool_type), intent(in) :: forcingPool, diagnosticsPool
      character(len=512) :: msg

      integer :: k

      real(kind=RKIND) :: &
         mu,                           &
         little_h_edge,                & ! little_h (fast-filtered BLD) at edge - used in streamfunction strength
         big_H_edge,                   & ! big_H (additional slow-filtered of little_h, an "effective" MLD)
         surfaceBuoyancyForcingEdge,   & ! buoyancy flux interpolated to edge
         tauXEdge,                     & ! zonal wind stress at edge
         tauYEdge,                     & ! meridional wind stress at edge
         windStressMagEdge,            & ! |tau| at edge
         u_starEdge,                   & ! friction velocity
         turbVelScale2,                & ! turbulent velocity scale denominator (filtered, interpolated from cells)
         gradBuoyML_bigH,              & ! buoyancy gradient recomputed over big_H
         hPartial,                     & ! partial layer thickness for big_H integration
         hML_local,                    & ! local copy of hML for integration
         hAv                             ! layer thickness average for integration
      
      real (kind=RKIND), dimension(:), pointer :: &
         surfaceBuoyancyForcing,       & ! surface buoyancy flux 
         boundaryLayerDepth,           & ! boundary layer depth from KPP
         windStressZonal,              & ! zonal wind stress 
         windStressMeridional            ! meridional wind stress 

      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth',  boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceBuoyancyForcing', surfaceBuoyancyForcing)
      call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)

      ! Psi = [ (C_r * Delta_s * |f| * hbar * Hbar^2 * gradBuoyML) /
      !         (m_* u_*^3 + n_* w_*^3)^(2/3) ] * mu(z)

      ! -------------------------------------------------------------------------
      ! Mixed-layer time filtering (MOM6-style two-timescale running mean) as shown
      ! in their src/parameterizations/lateral/MOM_mixed_layer_restrat.F90
      !
      ! The boundary layer depth used in the Bodner et al. (2023) MLE scheme is
      ! not used directly from the instantaneous KPP output. Instead, MOM6 applies
      ! a *direction-dependent running mean* to remove high-frequency variability
      ! (especially diurnal cycling) while preserving physically meaningful
      ! entrainment events.
      !
      ! This filter is implemented via rmean2ts(signal, filtered, tau_growing, tau_decaying, dt):
      !
      !   if signal >= filtered (deepening / entrainment):
      !       filtered <- (dt*signal + tau_growing*filtered) / (dt + tau_growing)
      !
      !   if signal < filtered (shoaling / restratification):
      !       filtered <- (dt*signal + tau_decaying*filtered) / (dt + tau_decaying)
      !
      ! Key properties:
      !
      !   - This is NOT is a regime-dependent relaxation with hysteresis.
      !   - Growing MLD (entrainment) responds more rapidly (smaller tau_growing).
      !   - Shoaling MLD is strongly damped (larger tau_decaying) to suppress
      !     high-frequency collapse of the mixed layer.
      !   - If tau = 0, the filter reduces to instantaneous resetting
      !     (used in older OM4-style implementations).
      !
      ! Physical intent:
      !   - Preserve deepening events associated with mixing and convection.
      !   - Suppress rapid shoaling due to transient stratification or diurnal noise.
      !   - Provide a dynamically smooth "active MLD" (little_h) for MLE physics.
      !
      ! Secondary filtering:
      !   A second application of the same two-timescale filter produces "big_H",
      !   a slower-evolving mixed-layer depth used in the MLE streamfunction
      !   scaling. This introduces additional temporal persistence at mesoscale
      !   time scales.
      !
      ! -------------------------------------------------------------------------
      ! Read pre-computed filtered fields (computed once per timestep in
      ! ocn_submesoscale_update_filtered_fields, not in this edge loop).
      ! This ensures consistent filtering across all cells and avoids race
      ! conditions from multiple edges updating the same cell.
      !
      ! Following MOM6 convention:
      !   - little_h (BLD_filtered) = fast-filtered BLD, represents the dynamically
      !     active boundary layer depth, used in the streamfunction strength (h * H^2)
      !   - big_H (BLD_filtered_slow) = slow-filtered MLD, provides temporal
      !     persistence at mesoscale timescales, used for H^2 scaling and shape function
      ! -------------------------------------------------------------------------

      ! Interpolate little_h (fast-filtered BLD) from cells to edge
      ! This is used for the "h" term in the streamfunction numerator
      little_h_edge = 0.5_RKIND * (BLD_filtered(cell1) + BLD_filtered(cell2))

      ! Interpolate big_H (slow-filtered MLD) from cells to edge
      ! This is used for the "H^2" term and the shape function mu(z)
      if (MLD_growing_Tfilt > 0.0_RKIND .or. MLD_decaying_Tfilt > 0.0_RKIND) then
         big_H_edge = 0.5_RKIND * (BLD_filtered_slow(cell1) + BLD_filtered_slow(cell2))
      else
         ! No slow filtering: use little_h for big_H as well
         big_H_edge = little_h_edge
      end if

      ! Interpolate surface buoyancy flux from centers to edge
      surfaceBuoyancyForcingEdge = 0.5_RKIND * (surfaceBuoyancyForcing(cell1) + & 
                                                surfaceBuoyancyForcing(cell2))

      ! Interpolate wind stress from centers to edge, compute magnitude and friction velocity.
      tauXEdge = 0.5_RKIND * (windStressZonal(cell1) + windStressZonal(cell2))
      tauYEdge = 0.5_RKIND * (windStressMeridional(cell1) + windStressMeridional(cell2))
      windStressMagEdge = (tauXEdge**2 + tauYEdge**2)**(0.5_RKIND)

      ! Compute friction velocity and impose a numerical floor for weakly forced regimes.
      u_starEdge = (windStressMagEdge/rho_sw)**(0.5_RKIND)
      u_starEdge = max(u_starEdge, 1e-18_RKIND)

      ! Read pre-computed filtered turbulent velocity scale from cells to edge
      ! (computed once per timestep in ocn_submesoscale_update_filtered_fields)
      turbVelScale2 = 0.5_RKIND * (turbVelScale2_filtered(cell1) + turbVelScale2_filtered(cell2))
      if (turbVelScale2 /= turbVelScale2) then
         streamFunction(:) = 0.0_RKIND
         return
      end if
      turbVelScale2 = max(turbVelScale2, 1.0E-18_RKIND)

      ! -------------------------------------------------------------------------
      ! Recompute buoyancy gradient over big_H depth for consistent vertical scaling
      !
      ! In MOM6, the streamfunction uses the slow-filtered mixed layer depth (big_H)
      ! for the H^2 scaling and the shape function mu(z). The buoyancy gradient
      ! is also recomputed as the thickness-weighted average over big_H, not over
      ! the instantaneous MLD.
      !
      ! We recompute the vertical average of gradDensityEddy from the surface
      ! down to big_H depth. This ensures consistency between the depth scale
      ! used in H^2 and the depth over which the buoyancy gradient is averaged.
      ! -------------------------------------------------------------------------

      ! Recompute buoyancy gradient averaged over big_H depth
      gradBuoyML_bigH = 0.0_RKIND
      hML_local = 0.0_RKIND
      do k = minLevelEdgeTop(iEdge)+1, maxLevelEdgeTop(iEdge)

         ! layer thickness (edge-mean, consistent with vertical integral elsewhere)
         hAv = layerThickEdgeMean_arr(k,iEdge)

         ! stop if we exceed big_H
         if (hML_local + hAv > big_H_edge) then
            hAv = big_H_edge - hML_local
         end if

         if (hAv <= 0.0_RKIND) exit

         gradBuoyML_bigH = gradBuoyML_bigH + hAv * gradDensityEddy_arr(k,iEdge)
         hML_local = hML_local + hAv

         if (hML_local >= big_H_edge) exit

      end do

      ! normalize
      if (big_H_edge <= 1.0E-12_RKIND .or. &
          big_H_edge /= big_H_edge) then
         if (.not. bodner_nan_reported) then
            write(msg,'(A,I6,2(A,I8),4(A,E14.6))') &
               'Bodner bad big_H at iEdge=', iEdge, &
               ' cell1=', cell1, &
               ' cell2=', cell2, &
               ' little_h=', little_h_edge, &
               ' big_H=', big_H_edge, &
               ' turbVel2=', turbVelScale2, &
               ' gradBuoy=', gradBuoyML_bigH
            call mpas_log_write(msg)
            bodner_nan_reported = .true.
         end if
         streamFunction(:) = 0.0_RKIND
         return
      end if
      gradBuoyML_bigH = gradBuoyML_bigH / max(big_H_edge, 1.0E-20_RKIND)
      ! convert to buoyancy gradient
      gradBuoyML_bigH = gravity * gradBuoyML_bigH / rho_sw

      do k = minLevelEdgeTop(iEdge)+1,maxLevelEdgeTop(iEdge)
         mu = compute_mle_shape(zEdge(k), big_H_edge)
         streamFunction(k) = ( (Cr*ds*abs(fEdge(iEdge))*little_h_edge * &
                                big_H_edge**2 * gradBuoyML_bigH) / &
                                turbVelScale2 ) * mu
      end do

  end subroutine compute_streamfunction_Bod2023

  !***********************************************************************

   !***********************************************************************
   !
   !  pure function rmean2ts
   !
   !> \brief
   !>   Two-timescale running mean filter (MOM6-style hysteretic filter).
   !>
   !> \details
   !>   If signal >= filtered: use tau_growing (allows rapid deepening).
   !>   If signal <  filtered: use tau_decaying (suppresses rapid shoaling).
   !>
   !>   Update: filtered <- (dt * signal + tau * filtered) / (dt + tau)
   !>
   !>   This creates direction-dependent relaxation with hysteresis,
   !>   useful for filtering high-frequency variability (e.g., diurnal
   !>   cycle) while preserving physical entrainment events.
   !>
   !>   Physical interpretation:
   !>     - tau_growing is the response time for a real deepening event.
   !>       A smaller tau_growing makes the filter follow entrainment more
   !>       quickly.
   !>     - tau_decaying is the smoothing time for shoaling or collapse of
   !>       the signal. A larger tau_decaying suppresses short-lived drops.
   !>
   !>   Example using the Bodner settings in this case:
   !>     BLD filter:
   !>       tau_growing  = 0 s
   !>       tau_decaying = 86400 s
   !>
   !>       - If the diagnosed boundary layer depth deepens from 20 m to 60 m,
   !>         the filtered BLD follows immediately because the growing time
   !>         scale is zero.
   !>       - If it then shoals back toward 20 m, the filter relaxes on a
   !>         1-day timescale, so short-lived shoaling pulses are damped.
   !>
   !>     MLD filter:
   !>       tau_growing  = 0 s
   !>       tau_decaying = 2.592e6 s
   !>
   !>       - A deepening MLD is again followed immediately.
   !>       - A brief shoaling pulse is retained much longer, with the filtered
   !>         MLD relaxing on a 30-day timescale.
   !>
   !>   Note: if tau_growing = tau_decaying = 0, the filter reduces to
   !>   instantaneous reset (filtered = signal).
   !
   !> \param[in] signal        - Unfiltered signal (arbitrary units)
   !> \param[in] filtered      - Current value of running mean (same units as signal)
   !> \param[in] tau_growing   - Time scale for growing signal [T ~> s]
   !> \param[in] tau_decaying  - Time scale for decaying signal [T ~> s]
   !> \param[in] dt            - Time step [T ~> s]
   !
   !***********************************************************************

   pure function rmean2ts(signal, filtered, tau_growing, tau_decaying, dt) &
       result(out)

       real(kind=RKIND), intent(in) :: signal
       real(kind=RKIND), intent(in) :: filtered
       real(kind=RKIND), intent(in) :: tau_growing
       real(kind=RKIND), intent(in) :: tau_decaying
       real(kind=RKIND), intent(in) :: dt

       real(kind=RKIND) :: out
       real(kind=RKIND) :: tau_eff, rt, aFac

       if (signal >= filtered) then
          tau_eff = tau_growing
       else
          tau_eff = tau_decaying
       end if

       rt = 1.0_RKIND / (dt + tau_eff)
       aFac = tau_eff * rt

       out = aFac * filtered + (1.0_RKIND - aFac) * signal

   end function rmean2ts

   !***********************************************************************
   !
   !  pure function compute_mle_shape
   !
   !> \brief
   !>   Computes the shape function mu(z) for the MLE streamfunction.
   !>
   !> \details
   !>   The shape function smoothly tapers the overturning circulation
   !>   to zero beneath the mixed layer. Uses the second-order shape
   !>   from Fox-Kemper et al. (2008, 2011):
   !>
   !>      mu = (1 - xi^2) * (1 + 5/21 * xi^2)
   !>
   !>   where xi = 2*z/zMLD + 1, with z=0 at the surface and z=-zMLD
   !>   at the base of the mixed layer.
   !
   !> \param[in] z     - Depth (positive downward, z <= 0) [Z ~> m]
   !> \param[in] zMLD  - Mixed layer depth [Z ~> m]
   !
   !***********************************************************************

   pure function compute_mle_shape(z, zMLD) result(mu)

       real(kind=RKIND), intent(in) :: z, zMLD

       real(kind=RKIND) :: mu, xi

       xi = 2.0_RKIND*z/zMLD + 1.0_RKIND
       mu = max(0.0_RKIND, (1.0_RKIND - xi**2) * (1.0_RKIND + 5.0_RKIND/21.0_RKIND * xi**2))

   end function compute_mle_shape

  !***********************************************************************
  !
  !  routine ocn_add_MLE_to_transport_vel
  !
  !> \brief   Submesoscale parameterization add to transport Vel
  !> \details
  !>  adds the MLE induced velocity to the transport velocity 
  !>  
  !***********************************************************************

  subroutine ocn_MLE_add_to_transport_vel(normalTransportVelocity, nEdges)

     real(kind=RKIND), dimension(:,:), intent(inout) :: &
         normalTransportVelocity

     integer, intent(in) :: nEdges
     
     integer :: iEdge, k

     if (.not. config_submesoscale_enable) return

     !$omp parallel
     !$omp do schedule(runtime) &
     !$omp private(k)
     do iEdge = 1, nEdges
        do k = 1, nVertLevels
           normalTransportVelocity(k,iEdge) = normalTransportVelocity(k,iEdge) + &
               normalMLEvelocity(k,iEdge)
        end do
     end do
     !$omp end do
     !$omp end parallel

  end subroutine ocn_MLE_add_to_transport_vel

  !***********************************************************************
  !
  !  routine ocn_submesoscale_init
  !
  !> \brief   Submesoscale parameterization init
  !> \details
  !>  Initializes parameters related to the submesoscale eddy parameterization 
  !>  
  !***********************************************************************

  subroutine ocn_submesoscale_init(err)

     integer, intent(out) :: err !< Output: error flag

     integer :: iEdge

     err = 0

     Ce = config_submesoscale_ce
     Lfmin = config_submesoscale_lfmin
     dsMax = config_submesoscale_ds_max
     tau = config_submesoscale_tau
     n_star = config_submesoscale_n_star
     m_star = config_submesoscale_m_star
     Cr = config_submesoscale_cr
     BLD_growing_Tfilt = config_submesoscale_bld_growing_tfilt
     BLD_decaying_Tfilt = config_submesoscale_bld_decaying_tfilt
     MLD_growing_Tfilt = config_submesoscale_mld_growing_tfilt
     MLD_decaying_Tfilt = config_submesoscale_mld_decaying_tfilt

     if (trim(config_submesoscale_scheme) /= 'FK2011' .and. &
         trim(config_submesoscale_scheme) /= 'Bod2023') then

         call mpas_log_write( &
            'Invalid value for config_submesoscale_scheme: ' // &
            trim(config_submesoscale_scheme), MPAS_LOG_CRIT)

         err = 1
         return

     end if

     allocate(time_scale(nEdgesAll))

     !$omp parallel
     !$omp do schedule(runtime)
     do iEdge=1,nEdgesAll
        time_scale(iEdge) = sqrt(fEdge(iEdge)**2 + tau**(-2.0_RKIND))
     end do
     !$omp end do
     !$omp end parallel

     ! Allocate storage for time-filtered BLD and MLD fields (only needed for Bodner scheme)
     if (trim(config_submesoscale_scheme) == 'Bod2023') then
         allocate(BLD_filtered(nCellsAll))
         allocate(BLD_filtered_slow(nCellsAll))
         allocate(turbVelScale2_filtered(nCellsAll))

         BLD_filtered(:) = 0.0_RKIND
         BLD_filtered_slow(:) = 0.0_RKIND
         turbVelScale2_filtered(:) = 0.0_RKIND

         filters_initialized = .false.
         bodner_nan_reported = .false.
      end if

   end subroutine ocn_submesoscale_init

   !***********************************************************************

   subroutine ocn_submesoscale_initialize_filters(diagnosticsPool, forcingPool)

      type(mpas_pool_type), intent(in) :: diagnosticsPool
      type(mpas_pool_type), intent(in) :: forcingPool

      integer :: iCell

      real(kind=RKIND) :: &
         boundaryLayerDepthCell, &
         surfaceBuoyancyForcingCell, &
         windStressMagCell, &
         u_starCell, &
         w_starCell_inst, &
         turbVelScale2_inst_cell

      real(kind=RKIND), dimension(:), pointer :: &
         boundaryLayerDepth, &
         surfaceBuoyancyForcing, &
         windStressZonal, &
         windStressMeridional

      call mpas_pool_get_array(diagnosticsPool, &
         'boundaryLayerDepth', boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, &
         'surfaceBuoyancyForcing', surfaceBuoyancyForcing)
      call mpas_pool_get_array(forcingPool, &
         'windStressZonal', windStressZonal)
      call mpas_pool_get_array(forcingPool, &
         'windStressMeridional', windStressMeridional)

      !$omp parallel do schedule(runtime) &
      !$omp private(iCell,boundaryLayerDepthCell,surfaceBuoyancyForcingCell, &
      !$omp         windStressMagCell,u_starCell,w_starCell_inst, &
      !$omp         turbVelScale2_inst_cell)
      do iCell = 1, nCellsAll
         boundaryLayerDepthCell = boundaryLayerDepth(iCell)

         ! Seed filters from instantaneous state
         BLD_filtered(iCell)      = boundaryLayerDepthCell
         BLD_filtered_slow(iCell) = boundaryLayerDepthCell

         surfaceBuoyancyForcingCell = surfaceBuoyancyForcing(iCell)

         windStressMagCell = sqrt( &
            windStressZonal(iCell)**2 + &
            windStressMeridional(iCell)**2 )

         u_starCell = sqrt(max(windStressMagCell/rho_sw,0.0_RKIND))

         w_starCell_inst = ( &
            max(0.0_RKIND,-surfaceBuoyancyForcingCell) * &
            boundaryLayerDepthCell )**(1.0_RKIND/3.0_RKIND)

         turbVelScale2_inst_cell = &
            (m_star*u_starCell**3 + &
               n_star*w_starCell_inst**3)**(2.0_RKIND/3.0_RKIND)

         turbVelScale2_filtered(iCell) = &
            max(turbVelScale2_inst_cell,1.0E-9_RKIND)

      end do
      !$omp end parallel do

   end subroutine ocn_submesoscale_initialize_filters

   !***********************************************************************
   !
   !  routine ocn_submesoscale_update_filtered_fields
   !
   !> \brief   Update time-filtered fields for Bodner et al. (2023) parameterization
   !> \details
   !>  This routine updates the filtered boundary layer depth and turbulent velocity
   !>  scale fields. It is called ONCE per timestep before the edge loop to ensure
   !>  consistent filtering across all cells.
   !>  
   !***********************************************************************

   subroutine ocn_submesoscale_update_filtered_fields(diagnosticsPool, forcingPool, dt)

      type(mpas_pool_type), intent(in) :: diagnosticsPool, forcingPool
      real(kind=RKIND), intent(in) :: dt

      integer :: iCell
      real(kind=RKIND) :: u_starCell, w_starCell_inst, turbVelScale2_inst_cell
      real(kind=RKIND) :: surfaceBuoyancyForcingCell, boundaryLayerDepthCell
      character(len=512) :: msg
      real(kind=RKIND) :: windStressZonalCell, windStressMeridionalCell, windStressMagCell
      real(kind=RKIND), dimension(:), pointer :: &
         boundaryLayerDepth,   &
         surfaceBuoyancyForcing, &
         windStressZonal,      &
         windStressMeridional

      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceBuoyancyForcing', surfaceBuoyancyForcing)
      call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)

      !------------------------------------------------------------------
      ! First-call initialization.
      !
      ! Seed the filters from the instantaneous state instead of zero.
      ! This avoids a completely artificial spinup of little_h and big_H.
      !------------------------------------------------------------------

      if (.not. filters_initialized) then

         !$omp parallel do schedule(runtime) &
         !$omp private(iCell, u_starCell, w_starCell_inst, turbVelScale2_inst_cell, &
         !$omp         surfaceBuoyancyForcingCell, boundaryLayerDepthCell, &
         !$omp         windStressZonalCell, windStressMeridionalCell, windStressMagCell)
         do iCell = 1, nCellsAll

            boundaryLayerDepthCell = boundaryLayerDepth(iCell)
            surfaceBuoyancyForcingCell = surfaceBuoyancyForcing(iCell)

            windStressZonalCell = windStressZonal(iCell)
            windStressMeridionalCell = windStressMeridional(iCell)

            windStressMagCell = sqrt( &
               windStressZonalCell**2 + &
               windStressMeridionalCell**2 )

            u_starCell = sqrt(windStressMagCell/rho_sw)
            u_starCell = max(u_starCell, 1.0e-18_RKIND)

            w_starCell_inst = ( &
               max(0.0_RKIND,-surfaceBuoyancyForcingCell) * &
               boundaryLayerDepthCell )**(1.0_RKIND/3.0_RKIND)

            turbVelScale2_inst_cell = &
               (m_star*u_starCell**3 + &
                  n_star*w_starCell_inst**3)**(2.0_RKIND/3.0_RKIND)

            turbVelScale2_inst_cell = &
               max(turbVelScale2_inst_cell, 1.0e-9_RKIND)

            ! Seed filters directly from instantaneous values

            BLD_filtered(iCell)      = boundaryLayerDepthCell
            BLD_filtered_slow(iCell) = boundaryLayerDepthCell

            turbVelScale2_filtered(iCell) = turbVelScale2_inst_cell

         end do
         !$omp end parallel do

         filters_initialized = .true.

         call mpas_log_write( &
            'Initialized Bod2023 filtered fields from instantaneous state.')

         return

      end if

      !$omp parallel do schedule(runtime) &
      !$omp private(iCell, u_starCell, w_starCell_inst, turbVelScale2_inst_cell, &
      !$omp         surfaceBuoyancyForcingCell, boundaryLayerDepthCell, &
      !$omp         windStressZonalCell, windStressMeridionalCell, windStressMagCell)
      do iCell = 1, nCellsAll
         boundaryLayerDepthCell = boundaryLayerDepth(iCell)
         surfaceBuoyancyForcingCell = surfaceBuoyancyForcing(iCell)

         ! Stage 1: Fast filter on instantaneous BLD -> little_h (BLD_filtered)
         if (BLD_growing_Tfilt > 0.0_RKIND .or. BLD_decaying_Tfilt > 0.0_RKIND) then
            BLD_filtered(iCell) = rmean2ts(boundaryLayerDepthCell, BLD_filtered(iCell), &
                                             BLD_growing_Tfilt, BLD_decaying_Tfilt, dt)
         else
            BLD_filtered(iCell) = boundaryLayerDepthCell
         end if

         ! Stage 2: Slow filter on little_h -> big_H (BLD_filtered_slow)
         if (MLD_growing_Tfilt > 0.0_RKIND .or. MLD_decaying_Tfilt > 0.0_RKIND) then
            BLD_filtered_slow(iCell) = rmean2ts(BLD_filtered(iCell), BLD_filtered_slow(iCell), &
                                                  MLD_growing_Tfilt, MLD_decaying_Tfilt, dt)
         else
            BLD_filtered_slow(iCell) = BLD_filtered(iCell)
         end if

         ! Compute instantaneous turbulent velocity scale for this cell
         windStressZonalCell = windStressZonal(iCell)
         windStressMeridionalCell = windStressMeridional(iCell)
         windStressMagCell = (windStressZonalCell**2 + windStressMeridionalCell**2)**(0.5_RKIND)
         u_starCell = (windStressMagCell/rho_sw)**(0.5_RKIND)
         u_starCell = max(u_starCell, 1e-18_RKIND)

         w_starCell_inst = (max(0.0_RKIND, -surfaceBuoyancyForcingCell) * &
                            boundaryLayerDepthCell)**(1.0_RKIND/3.0_RKIND)

         turbVelScale2_inst_cell = (m_star*u_starCell**3 + n_star*w_starCell_inst**3)**(2.0_RKIND/3.0_RKIND)
         turbVelScale2_inst_cell = max(turbVelScale2_inst_cell, 1e-9_RKIND)

         ! Filter turbulent velocity scale
         if (BLD_growing_Tfilt > 0.0_RKIND .or. BLD_decaying_Tfilt > 0.0_RKIND) then
            turbVelScale2_filtered(iCell) = rmean2ts(turbVelScale2_inst_cell, turbVelScale2_filtered(iCell), &
                                                       BLD_growing_Tfilt, BLD_decaying_Tfilt, dt)
         else
            turbVelScale2_filtered(iCell) = turbVelScale2_inst_cell
         end if

      end do
      !$omp end parallel do

   end subroutine ocn_submesoscale_update_filtered_fields

   !***********************************************************************

   subroutine ocn_submesoscale_finalize()

     deallocate(time_scale)
     
     if (allocated(BLD_filtered)) deallocate(BLD_filtered)
     if (allocated(BLD_filtered_slow)) deallocate(BLD_filtered_slow)
     if (allocated(turbVelScale2_filtered)) deallocate(turbVelScale2_filtered)
     bodner_nan_reported = .false.

   end subroutine ocn_submesoscale_finalize

end module ocn_submesoscale_eddies
