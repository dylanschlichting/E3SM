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
   real(kind=RKIND),dimension(:),allocatable :: time_scale
  !***********************************************************************

  contains

  !***********************************************************************
  !
  !  routine ocn_submesoscale_compute_velocity
  !
  !> \brief   Computes velocity from submesoscale eddies
  !> \details
  !>  This routine implements the Fox-Kemper et al 2011 submesoscale 
  !>  parameterization (https://doi.org/10.1016/j.ocemod.2010.09.002).
  !>  The transport velocity from the submesoscales is given by a 
  !>  stream function which is given as
  !>
  !>  Psi_FK = C_e*Delta_S/L_f*H^2 int(bGrad_H)/sqrt(f^2+tau^-2)*mu(z)
  !>
  !>  here C_e is a specified efficiency, Delta_s is taken as 
  !>  min(dcEdge,dsMax), H is the mixed layer depth, L_f is the subgrid
  !>  frontal width that is scaled to the coarse grid, f is the 
  !>  coriolis parameter, and tau is a time scale specified to prevent 
  !>  the singularity near the equator.  the horizontal buoyancy gradient
  !>  is integrated over the mixed layer depth and mu(z) is a non-dimensional
  !>  shape function that is zero below the mixed layer depth. 
  !> 
  !>  Psi_Bod = [ (C_r * Delta_s * |f| * h * H^2 * gradBuoyML) / 
  !>              (m_* u_*^3 + n_* w_*^3)^(2/3) ] * mu(z)
  !>
  !>  where C_r = C_e / C_L with C_L = (c_2)^2 * c_3 * Ri = 0.25 (Sec. 3),
  !>  h is the boundary layer depth (from KPP), u_* = (WindStressMag/rho_0)**(1/2) 
  !>  is the friction
  !>  velocity, and w_* is the convective velocity scale. This replaces
  !>  the old Lf (based on N, H, f) with a turbulence-based arrested 
  !>  frontal width that widens under strong surface forcing.
  !>  Note: u_* = (WindStressMag/rho_0)**(1/2)

  !***********************************************************************

  subroutine ocn_submesoscale_compute_velocity(forcingPool,diagnosticsPool)

     !-----------------------------------------------------------------
     !
     ! local variables
     !
     !-----------------------------------------------------------------

     type(mpas_pool_type), intent(in) :: forcingPool, diagnosticsPool
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
                diagnosticsPool)
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
      iEdge, cell1, cell2, zMLD, ds, gradBuoyML, zEdge, streamFunction, forcingPool, diagnosticsPool)

      integer, intent(in) :: iEdge, cell1, cell2
      real(kind=RKIND), intent(in) :: zMLD, ds, gradBuoyML
      real(kind=RKIND), dimension(:), intent(in) :: zEdge
      real(kind=RKIND), dimension(:), intent(inout) :: streamFunction
      type(mpas_pool_type), intent(in) :: forcingPool, diagnosticsPool

      integer :: k

      real(kind=RKIND) :: &
         mu,                           &
         Cr,                           & ! C_r = C_e / C_L
         boundaryLayerDepthEdge,       & ! boundary layer depth at edge
         surfaceBuoyancyForcingEdge,   & ! buoyancy flux interpolated to edge
         tauXEdge,                     & ! zonal wind stress at edge
         tauYEdge,                     & ! meridional wind stress at edge
         windStressMagEdge,            & ! |tau| at edge
         u_starEdge,                   & ! friction velocity
         w_starEdge,                   & ! convective velocity scale
         vel3,                         & ! m_* u_*^3 + n_* w_*^3
         turbVelScale                    ! turbulent velocity scale denominator
      
      real (kind=RKIND), dimension(:), pointer :: &
         surfaceBuoyancyForcing, & ! surface buoyancy flux 
         boundaryLayerDepth, &     ! boundary layer depth from KPP
         windStressZonal, &        ! zonal wind stress 
         windStressMeridional      ! meridional wind stress 

      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth',  boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceBuoyancyForcing', surfaceBuoyancyForcing)
      call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)

      ! Psi = [ (C_r * Delta_s * |f| * h * H^2 * gradBuoyML) /
      !         (m_* u_*^3 + n_* w_*^3)^(2/3) ] * mu(z)

      ! Assume 
      Cr = Ce / 0.25_RKIND
      ! Interpolate surface buoyancy flux and boundary layer depth from centers to edge
      surfaceBuoyancyForcingEdge = 0.5_RKIND * (surfaceBuoyancyForcing(cell1) + surfaceBuoyancyForcing(cell2))
      boundaryLayerDepthEdge = 0.5_RKIND * (boundaryLayerDepth(cell1) + boundaryLayerDepth(cell2))

      ! From registry, positive flux increases buoyancy, so use max(0, -flux) 
      ! to compute convective velocity scale
      w_starEdge = (max(0.0_RKIND, -surfaceBuoyancyForcingEdge) * boundaryLayerDepthEdge)**(1.0_RKIND/3.0_RKIND)

      ! Interpolate wind stress from centers to edge, compute magnitude and friction velocity
      tauXEdge = 0.5_RKIND * (windStressZonal(cell1) + windStressZonal(cell2))
      tauYEdge = 0.5_RKIND * (windStressMeridional(cell1) + windStressMeridional(cell2))
      windStressMagEdge = (tauXEdge**2 + tauYEdge**2)**(0.5_RKIND)

      ! Compute friction velocity 
      u_starEdge = (windStressMagEdge/rho_sw)**(0.5_RKIND)
      turbVelScale = (m_star*u_starEdge**3 + n_star*w_starEdge**3)**(2.0_RKIND/3.0_RKIND)
      ! To prevent a divide by zero error during very weak forcing, set a minimum value for the turbulent velocity
      ! scale is. 
      turbVelScale = max(turbVelScale, 1e-6_RKIND) 

      do k = minLevelEdgeTop(iEdge)+1,maxLevelEdgeTop(iEdge)
         mu = compute_mle_shape(zEdge(k), zMLD)
         streamFunction(k) = ( (Cr*ds*abs(fEdge(iEdge))*boundaryLayerDepthEdge* &
            zMLD**2*gradBuoyML) / turbVelScale ) * mu

      end do

  end subroutine compute_streamfunction_Bod2023

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

   end subroutine ocn_submesoscale_init

   !***********************************************************************

   subroutine ocn_submesoscale_finalize()

     deallocate(time_scale)

   end subroutine ocn_submesoscale_finalize

end module ocn_submesoscale_eddies