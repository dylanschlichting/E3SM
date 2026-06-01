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
             ocn_MLE_add_to_transport_vel

   real(kind=RKIND) :: tau, LfMin, Ce, dsMax, m_star, n_star
   real(kind=RKIND) :: BLD_growing_Tfilt, BLD_decaying_Tfilt
   real(kind=RKIND) :: MLD_growing_Tfilt, MLD_decaying_Tfilt
   real(kind=RKIND),dimension(:),allocatable :: time_scale
   
   ! Storage for time-filtered boundary layer and mixed layer depths
   ! These are updated each time step to implement direction-dependent running mean
   real(kind=RKIND),dimension(:),allocatable :: BLD_filtered      ! Running mean of BLD (fast)
   real(kind=RKIND),dimension(:),allocatable :: BLD_filtered_slow  ! Running mean of BLD (slow)
   real(kind=RKIND),dimension(:),allocatable :: wpup_filtered      ! Running mean of turbulent velocity scale
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
  !>   both wind-driven and convective turbulent forcing.
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

      ! Shared / FK2011 variables
     real(kind=RKIND) :: zMLD, bvfML, hML, gradBuoyML, bvfAv, hAv, ds

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
                iEdge,          &
                cell1,          &
                cell2,          &
                zMLD,           &
                ds,             &
                gradBuoyML,     &
                zEdge,          &
                streamFunction, & 
                forcingPool,    &
                diagnosticsPool, &
                dt)
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
      iEdge, cell1, cell2, zMLD, ds, gradBuoyML, zEdge, streamFunction, &
      forcingPool, diagnosticsPool, dt)

      integer, intent(in) :: iEdge, cell1, cell2
      real(kind=RKIND), intent(in) :: zMLD, ds, gradBuoyML, dt
      real(kind=RKIND), dimension(:), intent(in) :: zEdge
      real(kind=RKIND), dimension(:), intent(inout) :: streamFunction
      type(mpas_pool_type), intent(in) :: forcingPool, diagnosticsPool

      integer :: k

      real(kind=RKIND) :: &
         mu,                           &
         Cr,                           & ! C_r = C_e / C_L
         boundaryLayerDepthEdge,       & ! boundary layer depth at edge
         boundaryLayerDepthEdge_inst,  & ! instantaneous (unfiltered) BLD at edge
         surfaceBuoyancyForcingEdge,   & ! buoyancy flux interpolated to edge
         tauXEdge,                     & ! zonal wind stress at edge
         tauYEdge,                     & ! meridional wind stress at edge
         windStressMagEdge,            & ! |tau| at edge
         u_starEdge,                   & ! friction velocity
         w_starEdge,                   & ! convective velocity scale (fast-filtered BLD)
         w_starEdge_inst,              & ! convective velocity scale (instantaneous BLD)
         vel3,                         & ! m_* u_*^3 + n_* w_*^3
         turbVelScale2,                & ! turbulent velocity scale denominator (filtered)
         turbVelScale2_inst,           & ! turbulent velocity scale (instantaneous)
         BLD_input_cell1,              & ! BLD input for second-stage filter
         BLD_input_cell2                 ! BLD input for second-stage filter
      
      real (kind=RKIND), dimension(:), pointer :: &
         surfaceBuoyancyForcing,       & ! surface buoyancy flux 
         boundaryLayerDepth,           & ! boundary layer depth from KPP
         windStressZonal,              & ! zonal wind stress 
         windStressMeridional            ! meridional wind stress 

      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth',  boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceBuoyancyForcing', surfaceBuoyancyForcing)
      call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)

      ! Psi = [ (C_r * Delta_s * |f| * h * H^2 * gradBuoyML) /
      !         (m_* u_*^3 + n_* w_*^3)^(2/3) ] * mu(z)

      ! -------------------------------------------------------------------------
      ! Mixed-layer time filtering (MOM6-style two-timescale running mean) as shown
      ! in src/parameterizations/lateral/MOM_mixed_layer_restrat.F90
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
      !   - This is NOT a single exponential filter; it is a regime-dependent
      !     relaxation with hysteresis.
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

      ! Cr is the rescaled efficiency coefficient (C_e / C_L where C_L = 0.25)
      Cr = Ce / 0.25_RKIND

      ! Interpolate surface buoyancy flux from centers to edge
      surfaceBuoyancyForcingEdge = 0.5_RKIND * (surfaceBuoyancyForcing(cell1) + & 
                                                surfaceBuoyancyForcing(cell2))

      ! -------------------------------------------------------------------------
      ! MOM6-style two-stage filtering for Bodner et al. (2023) parameterization
      !
      ! Stage 1 (fast filter): Apply two-timescale running mean to instantaneous BLD
      !   to obtain "little_h" - the active mixing layer depth with diurnal cycle removed.
      !   Uses BLD_growing/decaying timescales.
      !
      ! Stage 2 (slow filter): Apply two-timescale running mean to little_h
      !   to obtain "big_H" - the mixed layer depth for streamfunction scaling.
      !   Uses MLD_growing/decaying timescales (typically slower).
      !
      ! This follows MOM6 convention where:
      !   - little_h represents the dynamically active boundary layer depth
      !   - big_H provides temporal persistence at mesoscale timescales
      ! -------------------------------------------------------------------------

      ! Stage 1: Fast filter on instantaneous BLD -> little_h
      if (BLD_growing_Tfilt > 0.0_RKIND .or. BLD_decaying_Tfilt > 0.0_RKIND) then
         ! Apply two-timescale filter at cell centers
         BLD_filtered(cell1) = rmean2ts(boundaryLayerDepth(cell1), BLD_filtered(cell1), &
                                          BLD_growing_Tfilt, BLD_decaying_Tfilt, dt)
         BLD_filtered(cell2) = rmean2ts(boundaryLayerDepth(cell2), BLD_filtered(cell2), &
                                          BLD_growing_Tfilt, BLD_decaying_Tfilt, dt)
      else
         ! No filtering: little_h = instantaneous BLD
         BLD_filtered(cell1) = boundaryLayerDepth(cell1)
         BLD_filtered(cell2) = boundaryLayerDepth(cell2)
      end if

      ! Stage 2: Slow filter on little_h -> big_H (if configured)
      if (MLD_growing_Tfilt > 0.0_RKIND .or. MLD_decaying_Tfilt > 0.0_RKIND) then
         ! Apply two-timescale filter to the already-filtered little_h
         BLD_filtered_slow(cell1) = rmean2ts(BLD_filtered(cell1), BLD_filtered_slow(cell1), &
                                                MLD_growing_Tfilt, MLD_decaying_Tfilt, dt)
         BLD_filtered_slow(cell2) = rmean2ts(BLD_filtered(cell2), BLD_filtered_slow(cell2), &
                                                MLD_growing_Tfilt, MLD_decaying_Tfilt, dt)
         ! Use big_H (slow-filtered) for streamfunction scaling
         boundaryLayerDepthEdge = 0.5_RKIND * (BLD_filtered_slow(cell1) + BLD_filtered_slow(cell2))
      else
         ! No slow filtering: use little_h (fast-filtered) for streamfunction scaling
         boundaryLayerDepthEdge = 0.5_RKIND * (BLD_filtered(cell1) + BLD_filtered(cell2))
      end if

      ! Interpolate wind stress from centers to edge, compute magnitude and friction velocity.
      ! This must be done BEFORE computing turbVelScale2, as u_starEdge is used in that calculation.
      tauXEdge = 0.5_RKIND * (windStressZonal(cell1) + windStressZonal(cell2))
      tauYEdge = 0.5_RKIND * (windStressMeridional(cell1) + windStressMeridional(cell2))
      windStressMagEdge = (tauXEdge**2 + tauYEdge**2)**(0.5_RKIND)

      ! Compute friction velocity and impose a numerical floor for weakly forced regimes.
      ! In MOM6, ustar_min_dflt = 2.0e-4 * omega * (GV%Angstrom_Z + GV%dZ_subroundoff)
      ! where omega = 7.291e-5 s^-1, Angstrom_Z ≈ 1e-10 m, and dZ_subroundoff is an
      ! even smaller roundoff-scale quantity. This gives ustar_min ~ O(10^-18 m s^-1),
      ! effectively zero and intended only to prevent exact division-by-zero.
      u_starEdge = (windStressMagEdge/rho_sw)**(0.5_RKIND)
      u_starEdge = max(u_starEdge, 1e-18_RKIND)

      ! From registry, positive flux increases buoyancy, but convection is driven by buoyancy 
      ! loss, so we flip the sign and then clip to keep only destabilizing contributions.
      ! MOM6 does NOT actually use w_star^2 directly in its MLE formulation.
      ! Instead, it constructs a turbulent velocity scale via:
      !   wpup = (m_*u_*^3 + n_*w_*^3)^(2/3)
      ! and then applies a floor (MIN_WSTAR2) to that quantity.
      ! The cube-root form here (w_starEdge = (-B0*h)^(1/3)) is only an intermediate
      ! scaling used to build w_star^3, not a prognostic or limited w_star^2 variable.

      ! First, compute wpup from INSTANTANEOUS (unfiltered) BLD and buoyancy flux
      boundaryLayerDepthEdge_inst = 0.5_RKIND * (boundaryLayerDepth(cell1) + boundaryLayerDepth(cell2))
      w_starEdge_inst = (max(0.0_RKIND, -surfaceBuoyancyForcingEdge) * &
                         boundaryLayerDepthEdge_inst)**(1.0_RKIND/3.0_RKIND)
      
      ! Compute unfiltered turbVelScale2
      turbVelScale2_inst = (m_star*u_starEdge**3 + n_star*w_starEdge_inst**3)**(2.0_RKIND/3.0_RKIND)
      ! In MOM6 examples, w'u' is capped at 1e-9 to prevent numerical issues even though their 
      ! default values are 1e-24. 
      turbVelScale2_inst = max(turbVelScale2_inst, 1e-9_RKIND)
      
      ! Apply two-timescale filter to turbVelScale2 using fast timescales (BLD_growing/decaying)
      ! This follows MOM6 convention: wpup is filtered with the same fast timescales used for BLD
      if (BLD_growing_Tfilt > 0.0_RKIND .or. BLD_decaying_Tfilt > 0.0_RKIND) then
         ! Filter at cell centers and then average to edge
         wpup_filtered(cell1) = rmean2ts(turbVelScale2_inst, wpup_filtered(cell1), &
                                           BLD_growing_Tfilt, BLD_decaying_Tfilt, dt)
         wpup_filtered(cell2) = rmean2ts(turbVelScale2_inst, wpup_filtered(cell2), &
                                           BLD_growing_Tfilt, BLD_decaying_Tfilt, dt)
         turbVelScale2 = 0.5_RKIND * (wpup_filtered(cell1) + wpup_filtered(cell2))
      else
         turbVelScale2 = turbVelScale2_inst
      end if

      ! Now compute w_starEdge for use in the streamfunction (using filtered BLD for consistency)
      ! NOTE: Use boundaryLayerDepthEdge which is either the fast-filtered BLD (if
      ! MLD_growing/decaying_Tfilt are zero) or the slow-filtered BLD (if they are non-zero).
      w_starEdge = (max(0.0_RKIND, -surfaceBuoyancyForcingEdge) * & 
                        boundaryLayerDepthEdge)**(1.0_RKIND/3.0_RKIND)

      do k = minLevelEdgeTop(iEdge)+1,maxLevelEdgeTop(iEdge)
         mu = compute_mle_shape(zEdge(k), zMLD)
         streamFunction(k) = ( (Cr*ds*abs(fEdge(iEdge))*boundaryLayerDepthEdge* &
            zMLD**2*gradBuoyML) / turbVelScale2 ) * mu
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
        allocate(wpup_filtered(nCellsAll))

        BLD_filtered(:) = 0.0_RKIND
        BLD_filtered_slow(:) = 0.0_RKIND
        wpup_filtered(:) = 0.0_RKIND
     end if

   end subroutine ocn_submesoscale_init

   !***********************************************************************

   subroutine ocn_submesoscale_finalize()

     deallocate(time_scale)
     
     if (allocated(BLD_filtered)) deallocate(BLD_filtered)
     if (allocated(BLD_filtered_slow)) deallocate(BLD_filtered_slow)
     if (allocated(wpup_filtered)) deallocate(wpup_filtered)

   end subroutine ocn_submesoscale_finalize

end module ocn_submesoscale_eddies