! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_discrete_variance_decay
!
!> \brief MPAS ocean analysis mode member: Computes the physical and numerical 
!>        destruction of tracer variance for active tracers. 

!> \author Dylan Schlichting
!> \date   April-May 2025
!> \details
!>   Computes destruction of tracer variance via the numerical implementation
!>   of the 3D transport equation 
!>   Returns a nVertLevels x nCells x Time array for each tracer. 
!> 
!>   Details described in Burchard & Rennau (2008) Ocean Modelling
!>   https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JC008108
!>  
!>   Notes:
!>   chi_phy_vert = 2 * \kappa_v * (\partial_z c)**2
!>   chi_phy_horz = 2 * [ \kappa_h * (\partial_x c)**2 + (\partial_y c)**2 ] 
!>   chi_num_br08 = [ A{c**2} - (A{c})**2 ] / dt_baroclinic 
!>   chi_num_k14 = 
!>   Vertical gradients modified from Eliassen Palm analysis member,
!>   Vertical diffusion coefficient imported from vertical mixing 
!>   To do: Add support for explict horizontal mixing, which is not used 
!>   by default. Add GM + Redi support
!-----------------------------------------------------------------------

module ocn_discrete_variance_decay

    use mpas_derived_types
    use mpas_pool_routines
    use mpas_dmpar
    use ocn_config
    use mpas_timekeeping
    use mpas_stream_manager
 
    use ocn_constants
    use ocn_diagnostics_variables
    use ocn_vmix
    use ocn_tracer_advection_mono
 
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
 
    public :: ocn_init_discrete_variance_decay, &
              ocn_compute_discrete_variance_decay, &
              ocn_restart_discrete_variance_decay, &
              ocn_finalize_discrete_variance_decay, & 
              computeVerticalDerivative
 
    !--------------------------------------------------------------------
    !
    ! Private module variables
    !
    !--------------------------------------------------------------------
 
 !***********************************************************************
 
 contains
 
 !***********************************************************************
 !
 !  routine ocn_init_discrete_variance_decay
 !
 !> \brief   Initialize MPAS-Ocean analysis member
 !> \author Dylan Schlichting
 !> \date    Apr 07, 2025
 !> \details
 !>  This routine conducts all initializations required for the
 !>  MPAS-Ocean analysis member.
 !
 !-----------------------------------------------------------------------
 
    subroutine ocn_init_discrete_variance_decay(domain, err)!{{{
 
       !-----------------------------------------------------------------
       !
       ! input variables
       !
       !-----------------------------------------------------------------
 
       !-----------------------------------------------------------------
       !
       ! input/output variables
       !
       !-----------------------------------------------------------------
 
       type(domain_type), intent(inout) :: domain
 
       !-----------------------------------------------------------------
       !
       ! output variables
       !
       !-----------------------------------------------------------------
 
       integer, intent(out) :: err !< Output: error flag
 
       !-----------------------------------------------------------------
       !
       ! local variables
       !
       !-----------------------------------------------------------------
 
       err = 0
 
    end subroutine ocn_init_discrete_variance_decay!}}}
 
 !***********************************************************************
 !
 !  routine ocn_compute_discrete_variance_decay
 !
 !> \brief   Compute MPAS-Ocean analysis member
 !> \author  Dylan Schlichting
 !> \date    April-May 2025
 !> \details
 !>  This routine conducts all computation required for this
 !>  MPAS-Ocean analysis member.
 !
 !-----------------------------------------------------------------------
 
    subroutine ocn_compute_discrete_variance_decay(domain, timeLevel, err)!{{{
 
       !-----------------------------------------------------------------
       !
       ! input variables
       !
       !-----------------------------------------------------------------
 
       integer, intent(in) :: timeLevel
 
       !-----------------------------------------------------------------
       !
       ! input/output variables
       !
       !-----------------------------------------------------------------
 
       type(domain_type), intent(inout) :: domain
 
       !-----------------------------------------------------------------
       !
       ! output variables
       !
       !-----------------------------------------------------------------
 
       integer, intent(out) :: err !< Output: error flag
 
       !-----------------------------------------------------------------
       !
       ! local variables
       !
       !-----------------------------------------------------------------
 
       type(mpas_pool_type), pointer :: discreteVarianceDecayAMPool
       type(dm_info) :: dminfo
       type(block_type), pointer :: block
       type(mpas_pool_type), pointer :: statePool
       type(mpas_pool_type), pointer :: meshPool
       type(mpas_pool_type), pointer :: tracersPool
       type(mpas_pool_type), pointer :: diagnosticsPool

       real(kind=RKIND), dimension(:,:,:), pointer :: activeTracers
       real(kind=RKIND), dimension(:,:,:), pointer :: tracers
       real(kind=RKIND), dimension(:,:,:), allocatable :: activeTracersSquared
       real(kind=RKIND), dimension(:,:,:), allocatable :: squareAdvectedHorzTracer
       real(kind=RKIND), dimension(:,:,:), allocatable :: squareAdvectedVertTracer

       real(kind=RKIND), dimension(:,:,:), pointer :: chiSpurHorzTracerBR08, chiSpurVertTracerBR08

       ! pointers for tendencies used in diagnostic budget computation, A{c}
       real (kind=RKIND), dimension(:,:,:), pointer, contiguous :: &
             activeTracerHorizontalAdvectionTendency, &
             activeTracerVerticalAdvectionTendency,   &
             activeTracerHorizontalAdvectionTendencySquared, &
             activeTracerVerticalAdvectionTendencySquared
       
      real (kind=RKIND), dimension(:,:,:), pointer, contiguous :: &
             activeTracerSquaredHorizontalAdvectionTendency, &
             activeTracerSquaredVerticalAdvectionTendency
 
       ! Here are some example variables which may be needed for your analysis member
       real(kind=RKIND), dimension(:,:), pointer :: dTdz
       real(kind=RKIND), dimension(:,:), pointer :: dsdz
       real(kind=RKIND), dimension(:,:), pointer :: vertDiffTopOfCell
       real(kind=RKIND), dimension(:,:), pointer :: chiPhyVerSalt, chiPhyVerTemp
       real(kind=RKIND), dimension(:), pointer :: areaCell, bottomDepth
       real(kind=RKIND), pointer :: config_cvmix_background_diffusion
       real(kind=RKIND), allocatable :: vertDiffMid(:,:)
       real(kind=RKIND) :: backgroundDiff
       real(kind=RKIND), dimension(:,:), allocatable :: tmp2D
       real(kind=RKIND), dimension(:,:), allocatable :: sal2D
       
       integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot
       integer, pointer :: index_temperature, index_salinity, & 
       nVertLevels, nCells, nEdgesSolve, nVerticesSolve
       integer :: iTracer, k, iCell

       logical, pointer :: config_cvmix_background_diffusion_enable
       logical :: computeBudgets  = .false.

       character(len=StrKIND), pointer :: config_dt
       type (MPAS_timeInterval_type) :: timeStepESMF
       real (kind=RKIND) :: dt,

       err = 0
 
       dminfo = domain%dminfo

       call mpas_pool_get_config(configs, 'config_dt', config_dt)
       call mpas_set_timeInterval(timeStepESMF, timeString=config_dt,ierr=err)
       print *, 'dt = ', dt
       
 
       block => domain%blocklist
       do while (associated(block))
          
          call mpas_pool_get_subpool(block%structs, 'state', statePool)
          call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block%structs, 'diagnostics', diagnosticsPool)
          call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
          call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
          call mpas_pool_get_subpool(block%structs, 'discreteVarianceDecayAM', discreteVarianceDecayAMPool)

          call mpas_pool_get_dimension(block%dimensions, 'nVertLevels', nVertLevels)
          call mpas_pool_get_dimension(block%dimensions, 'nCells', nCells)

          call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
          call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
          call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
          call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
          call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
          call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
          call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

          call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
          call mpas_pool_get_array(discreteVarianceDecayAMPool, 'chiPhyVerSalt', chiPhyVerSalt) 
          call mpas_pool_get_array(discreteVarianceDecayAMPool, 'chiPhyVerTemp', chiPhyVerTemp) 
         !  if (config_use_br08_spur) then
         call mpas_pool_get_array(discreteVarianceDecayAMPool, 'chiSpurHorzTracerBR08', chiSpurHorzTracerBR08) 
         call mpas_pool_get_array(discreteVarianceDecayAMPool, 'chiSpurVertTracerBR08', chiSpurVertTracerBR08) 
         !  end if 
          !  if (config_use_k14)
         !  call mpas_pool_get_array(discreteVarianceDecayAMPool, 'chiSpurHorzTracerK14', chiSpurHorzTracerK14) 
         !  call mpas_pool_get_array(discreteVarianceDecayAMPool, 'chiSpurVertTracerK14', chiSpurVertTracerK14) 
          
          ! Variables needed for numerical mixing
          call mpas_pool_get_array(diagnosticsPool, 'activeTracerHorizontalAdvectionTendency', &
                                   activeTracerHorizontalAdvectionTendency)
          call mpas_pool_get_array(diagnosticsPool, 'activeTracerVerticalAdvectionTendency',  &
                                   activeTracerVerticalAdvectionTendency)

          ! Interpolate vertical diffusivity from top of cell to middle of layer
          allocate(vertDiffMid(nVertLevels-1, nCells))
          do iCell = 1, nCells
             do k = 1, nVertLevels-1
                 vertDiffMid(k,iCell) = ( vertDiffTopOfCell(k,iCell) + vertDiffTopOfCell(k+1,iCell)) / 2.0d0
             end do
          end do

          ! Computer vertical temperature gradient
          ! Initialize temp, salt, and their derivatives
          allocate(tmp2D(nVertLevels,nCells))
          allocate(dTdz(nVertLevels,nCells))
          tmp2D = activeTracers(index_temperature, :, :)
          call computeVerticalDerivative(nCells, nVertLevels, zMid, & 
                                         tmp2D, dTdz)

          allocate(sal2D(nVertLevels,nCells))
          allocate(dsdz(nVertLevels,nCells))
          sal2D = activeTracers(index_salinity, :, :)
          ! Computer vertical salinity gradient
          call computeVerticalDerivative(nCells, nVertLevels, zMid, & 
                                         sal2D, dsdz)

          ! Local estimates of vertical physical mixing 
          do iCell = 1, nCells
              do k = 1, maxLevelCell(iCell)
                  ! Subtract background tracer diffusion if turned on  
                  backgroundDiff = config_cvmix_background_diffusion 
                  ! Temperature 
                  chiPhyVerTemp(k,iCell) = 2.0_RKIND * (vertDiffMid(k,iCell)-backgroundDiff) & 
                  * ( (dTdz(k,iCell) * dTdz(k,iCell)) )
                  ! Salinity 
                  chiPhyVerSalt(k,iCell) = 2.0_RKIND * (vertDiffMid(k,iCell)-backgroundDiff) & 
                  * ( (dsdz(k,iCell) * dsdz(k,iCell)) )
              end do
          end do
          
          ! Local estimates of horizontal physical mixing due to laplacian diffusion 
          ! Local estimates of horizontal physical mixing due to biharmonic diffusion 
          ! Local estimates of redi diffusion 

          ! Local estimates of numerical mixing following Burchard & Rennau (2008) O.M.,
          ! (A{c^2} - (A{c})^2) / dt. Split into horizontal and vertical components, 
          ! Calculate (A{c})^2) first because it is already computed at run time.
         !  if (config_use_br08_spur) then 
         !  allocate(squareAdvectedHorzTracer(size(activeTracerHorizontalAdvectionTendency,1), size(activeTracerHorizontalAdvectionTendency,2), size(activeTracerHorizontalAdvectionTendency,3)))
         !  allocate(squareAdvectedVertTracer(size(activeTracerVerticalAdvectionTendency,1), size(activeTracerVerticalAdvectionTendency,2), size(activeTracerVerticalAdvectionTendency,3)))
          
         !  squareAdvectedHorzTracer(1,:,:) = activeTracerHorizontalAdvectionTendency(1,:,:)*activeTracerHorizontalAdvectionTendency(1,:,:)
         !  squareAdvectedVertTracer(1,:,:) = activeTracerVerticalAdvectionTendency(1,:,:)*activeTracerVerticalAdvectionTendency(1,:,:)

         ! !  print *, 'Surface horizontal advected temperature**2', activeTracerHorizontalAdvectionTendency(1,1,:)
         ! !  print *, 'Surface horizontal advected salinity**2', activeTracerHorizontalAdvectionTendency(2,1,:)

         !  squareAdvectedHorzTracer(2,:,:) = activeTracerHorizontalAdvectionTendency(2,:,:)*activeTracerHorizontalAdvectionTendency(2,:,:)
         !  squareAdvectedVertTracer(2,:,:) = activeTracerVerticalAdvectionTendency(2,:,:)*activeTracerVerticalAdvectionTendency(2,:,:)

         ! !  call mpas_timer_start("tracer adv")

         !  ! Square active tracers, then passive them through advection scheme. 
         !  ! call specific advection routine based on choice of monotonicity
         !  allocate(activeTracersSquared(size(activeTracers,1), size(activeTracers,2), size(activeTracers,3)))

         !  activeTracersSquared(1,:,:) = activeTracers(1,:,:) * activeTracers(1,:,:)
         !  activeTracersSquared(2,:,:) = activeTracers(2,:,:) * activeTracers(2,:,:)

         !  ! call specific advection routine based on choice of monotonicity, this updates activeTracerTendency components
         !  call ocn_tracer_advection_tend(activeTracersSquared, &
         !                                 normalThicknessFlux, vertAleTransportTop, &
         !                                 layerThickness, dt, tracerGroupTend, &
         !                                 computeBudgets)
         ! !  if (config_flux_limiter == 'monotonic') then
         ! !      call ocn_tracer_advection_mono_tend(tend, activeTracersSquared, layerThickness,    &
         ! !                                          normalThicknessFlux, w, dt, &
         ! !                                          computeBudgets)
         ! !  else
         ! !      call ocn_tracer_advection_std_tend(activeTracersSquared, normalThicknessFlux, &
         ! !                                         w, layerThickness, dt, tend)
         ! !  endif

         !  print *, 'Surface temperature**2 advected', activeTracerHorizontalAdvectionTendency(1,1,:) 
         !  print *, 'Surface salinity**2 advected', activeTracerHorizontalAdvectionTendency(2,1,:) 
      
         !  !  call mpas_timer_stop("tracer adv")
         !  ! Compute numerical mixing
         !  chiSpurHorzTracerBR08(1,:,:) = ( activeTracerHorizontalAdvectionTendency(1,:,:) - squareAdvectedHorizontalTracer(1,:,:) ) / dt
         !  chiSpurVertTracerBR08(1,:,:) = ( activeTracerVerticalAdvectionTendency(1,:,:) - squareAdvectedVerticalTracer(1,:,:) ) / dt

         !  chiSpurHorzTracerBR08(2,:,:) = ( activeTracerHorizontalAdvectionTendency(2,:,:) - squareAdvectedHorizontalTracer(2,:,:) ) / dt
         !  chiSpurVertTracerBR08(2,:,:) = ( activeTracerVerticalAdvectionTendency(2,:,:) - squareAdvectedVerticalTracer(2,:,:) ) / dt
         ! !  end if
          block => block%next
       end do
       deallocate(vertDiffMid)
       deallocate(tmp2D)
       deallocate(sal2D)
       deallocate(squareAdvectedHorzTracer)
       deallocate(squareAdvectedVertTracer)
       deallocate(activeTracersSquared)
 
    end subroutine ocn_compute_discrete_variance_decay!}}}

 !***********************************************************************
 !
 !  routine computeVerticalDerivative
 !
 !> \brief   Computes vertical gradient 
 !> \author  Dylan Schlichting
 !> \date    April-May 2025
 !> \details
 !>  Computes the vertical gradient of a scalar field.
 !>  Modified from mpas_ocn_eliassen_palm.F with help from ChatGPT
 !>  so that firstLayer is always 1 and the last layer is the total 
 !>  number of vertical layers.
 !-----------------------------------------------------------------------

    subroutine computeVerticalDerivative(nCells, nVertLevels, heightMid, field, derivativeField)
      integer, intent(in) :: nCells, nVertLevels
      real(kind=RKIND), dimension(:,:), intent(in) :: heightMid
      real(kind=RKIND), dimension(:,:), intent(in) :: field
      real(kind=RKIND), dimension(:,:), intent(out) :: derivativeField
      ! real(kind=RKIND), dimension(:,:), allocatable, intent(out) :: derivativeField

      integer :: iCell, kLayer
      integer :: firstLayer, lastLayer

      real(kind=RKIND) :: wrkAbove, wrkBelow, dz

      ! Set firstLayer to always be 1, and lastLayer to always be nVertLevels
      firstLayer = 1
      lastLayer = nVertLevels  ! Treat lastLayer as constant across all cells
  
      ! Initialize the derivativeField, note this is allocated before call time
      ! otherwise you well get a memory seg fault
      derivativeField = 0.0_RKIND
  
      do iCell = 1, nCells
          ! Calculate derivative at the first layer
          wrkAbove = field(firstLayer, iCell)
          wrkBelow = field(firstLayer + 1, iCell)
          dz = heightMid(firstLayer, iCell) - heightMid(firstLayer + 1, iCell)
          derivativeField(firstLayer, iCell) = (wrkAbove - wrkBelow) / dz  

          ! Loop over intermediate layers
          do kLayer = 2, lastLayer - 1
              wrkAbove = field(kLayer - 1, iCell)
              wrkBelow = field(kLayer + 1, iCell)
              dz = heightMid(kLayer - 1, iCell) - heightMid(kLayer + 1, iCell)
              derivativeField(kLayer, iCell) = (wrkAbove - wrkBelow) / dz
          end do
  
          ! Calculate derivative at the last layer
          wrkAbove = field(lastLayer - 1, iCell)
          wrkBelow = field(lastLayer, iCell)
          dz = heightMid(lastLayer - 1, iCell) - heightMid(lastLayer, iCell)
  
          derivativeField(lastLayer, iCell) = (wrkAbove - wrkBelow) / dz
      end do
  
  end subroutine computeVerticalDerivative
 
 !***********************************************************************
 !
 !  routine ocn_restart_discrete_variance_decay
 !
 !> \brief   Save restart for MPAS-Ocean analysis member
 !> \author  Dylan Schlichting
 !> \date    February 06, 2025
 !> \details
 !>  This routine conducts computation required to save a restart state
 !>  for the MPAS-Ocean analysis member.
 !
 !-----------------------------------------------------------------------
 
    subroutine ocn_restart_discrete_variance_decay(domain, err)!{{{
 
       !-----------------------------------------------------------------
       !
       ! input variables
       !
       !-----------------------------------------------------------------
 
       !-----------------------------------------------------------------
       !
       ! input/output variables
       !
       !-----------------------------------------------------------------
 
       type(domain_type), intent(inout) :: domain
 
       !-----------------------------------------------------------------
       !
       ! output variables
       !
       !-----------------------------------------------------------------
 
       integer, intent(out) :: err !< Output: error flag
 
       !-----------------------------------------------------------------
       !
       ! local variables
       !
       !-----------------------------------------------------------------
 
       err = 0
 
    end subroutine ocn_restart_discrete_variance_decay!}}}
 
 !***********************************************************************
 !
 !  routine ocn_finalize_discrete_variance_decay
 !
 !> \brief   Finalize MPAS-Ocean analysis member
 !> \author  Dylan Schlichting
 !> \date    Apr 07, 2025
 !> \details
 !>  This routine conducts all finalizations required for this
 !>  MPAS-Ocean analysis member.
 !
 !-----------------------------------------------------------------------
 
    subroutine ocn_finalize_discrete_variance_decay(domain, err)!{{{
 
       !-----------------------------------------------------------------
       !
       ! input variables
       !
       !-----------------------------------------------------------------
 
       !-----------------------------------------------------------------
       !
       ! input/output variables
       !
       !-----------------------------------------------------------------
 
       type(domain_type), intent(inout) :: domain
 
       !-----------------------------------------------------------------
       !
       ! output variables
       !
       !-----------------------------------------------------------------
 
       integer, intent(out) :: err !< Output: error flag
 
       !-----------------------------------------------------------------
       !
       ! local variables
       !
       !-----------------------------------------------------------------
 
       err = 0
 
    end subroutine ocn_finalize_discrete_variance_decay!}}}
 
end module ocn_discrete_variance_decay
 
 ! vim: foldmethod=marker