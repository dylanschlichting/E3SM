! Copyright (c) 2025,  Los Alamos National Security, LLC (LANS)
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
!>   Returns a nVertLevelsP1 x nCells x Time array for physical mixing.
!>   Returns a nVertLevels x nCells x Time array for numerical mixing. 
!> 
!>   Details described in Burchard & Rennau (2008) Ocean Modelling
!>   https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JC008108
!>  
!>   Notes:
!>   chi_phy_vert = 2 * \kappa_v * (\partial_z c)**2
!>   chi_phy_horz = 2 * [ \kappa_h * (\partial_x c)**2 + (\partial_y c)**2 ] 
!>   chi_num_br08 = [ A{c**2} - (A{c})**2 ] / dt_baroclinic 
!>   chi_num_k14 = 
!>   Vertical gradients using the approach in subroutine ocn_GM_compute_Bolus_velocity
!>   Vertical diffusion coefficient imported from vertical mixing 
!>   Numerical and physical located are not co-located on the vertical grid.
!>   Numerical mixing is on the mid-layer points, physical mixing is on cell tops
!>   because the vertical diffusivity is defined at the cell top.
!>   To do: Add support for explict horizontal mixing, which is not used 
!>   by default. Add GM + Redi support, much more complex. Add K14 method
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
    use mpas_log
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
              ocn_finalize_discrete_variance_decay
 
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
       real(kind=RKIND), dimension(:,:,:), pointer :: numericalMixingTracers
       real(kind=RKIND), dimension(:,:,:), pointer :: tracers
       real(kind=RKIND), dimension(:,:,:), pointer :: chiSpurTracerBR08

       ! Vert temp & salt grad, vert diffusivity
       real(kind=RKIND), dimension(:,:), pointer ::  vertDiffTopOfCell
       real(kind=RKIND), dimension(:,:), allocatable :: dTdzTop, dSdzTop, dzTop
       ! Vertical destruction of salt and temp variance
       real(kind=RKIND), dimension(:,:), pointer :: chiPhyVerSalt, chiPhyVerTemp
       real(kind=RKIND), dimension(:), pointer :: areaCell, bottomDepth

       real(kind=RKIND) :: backgroundDiff
       ! Temporary arrays for intermediate calcs
       real (kind=RKIND), dimension(:,:), pointer :: layerThickness
 
       integer, dimension(:), pointer :: minLevelCell, maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot
       integer, pointer :: index_temperature, index_salinity, & 
       nVertLevels, nVertLevelsP1, nCells, nCellsSolve, nEdgesSolve, nVerticesSolve
       integer :: iTracer, k, iCell

       character(len=StrKIND), pointer :: config_dt
       type (MPAS_timeInterval_type) :: timeStepESMF
       real (kind=RKIND) :: dt

       err = 0
 
       dminfo = domain%dminfo
       ! Convert dt to real
       call mpas_pool_get_config(domain%configs, 'config_dt', config_dt)
       call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err)
       call mpas_get_timeInterval(timeStepESMF, dt=dt)
 
       block => domain%blocklist
       do while (associated(block))
          call mpas_pool_get_subpool(block%structs, 'state', statePool)
          call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block%structs, 'diagnostics', diagnosticsPool)
          call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
          call mpas_pool_get_subpool(block%structs, 'discreteVarianceDecayAM', discreteVarianceDecayAMPool)

          call mpas_pool_get_dimension(block%dimensions, 'nVertLevels', nVertLevels)
          call mpas_pool_get_dimension(block%dimensions, 'nVertLevelsP1', nVertLevelsP1)
          call mpas_pool_get_dimension(block%dimensions, 'nCells', nCells)
          call mpas_pool_get_dimension(block%dimensions, 'nCellsSolve', nCellsSolve)

          call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
          call mpas_pool_get_array(tracersPool, 'numericalMixingTracers', numericalMixingTracers, 1)
          call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
          call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

          call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
          call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
          call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

          call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

          call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
          call mpas_pool_get_array(discreteVarianceDecayAMPool, 'chiSpurTracerBR08', chiSpurTracerBR08) 
          call mpas_pool_get_array(discreteVarianceDecayAMPool, 'chiPhyVerSalt', chiPhyVerSalt) 
          call mpas_pool_get_array(discreteVarianceDecayAMPool, 'chiPhyVerTemp', chiPhyVerTemp)

          ! Calculate vertical temperature and salinity gradients. Loop over layers 2 to n-1.
          ! Modified from components/mpas-ocean/src/shared/mpaso_ocn_gm.F
          allocate(dTdzTop(nVertLevels + 1, nCells))
          allocate(dSdzTop(nVertLevels + 1, nCells))
          allocate(dzTop(nVertLevels + 1, nCells))  
          do iCell = 1, nCellsSolve
            do k = minLevelCell(iCell)+1, maxLevelCell(iCell)
               dzTop(k,iCell) = 0.5_RKIND*(layerThickness(k - 1, iCell) + layerThickness(k, iCell))
               dTdzTop(k,iCell) = (activeTracers(index_temperature, k - 1, iCell) &
                                  - activeTracers(index_temperature, k, iCell)) &
                                  / dzTop(k,iCell)
               dSdzTop(k,iCell) = (activeTracers(index_salinity, k - 1, iCell) &
                                  - activeTracers(index_salinity, k, iCell)) &
                                  / dzTop(k,iCell)
            end do
            dzTop(1:minLevelCell(iCell),iCell) = 0.0_RKIND
            dTdzTop(1:minLevelCell(iCell),iCell) = 0.0_RKIND
            dSdzTop(1:minLevelCell(iCell),iCell) = 0.0_RKIND
            dzTop(maxLevelCell(iCell) + 1,iCell) = 0.0_RKIND
            dTdzTop(maxLevelCell(iCell) + 1,iCell) = 0.0_RKIND
            dSdzTop(maxLevelCell(iCell) + 1,iCell) = 0.0_RKIND
          end do

         ! Local estimates of vertical physical mixing. 
         do iCell = 1, nCellsSolve
              do k = minLevelCell(iCell), maxLevelCell(iCell)
                  ! Temperature 
                  chiPhyVerTemp(k,iCell) = 2.0_RKIND * (vertDiffTopOfCell(k,iCell)) & 
                  * ( (dTdzTop(k,iCell) * dTdzTop(k,iCell)) )
                  ! Salinity 
                  chiPhyVerSalt(k,iCell) = 2.0_RKIND * (vertDiffTopOfCell(k,iCell)) & 
                  * ( (dSdzTop(k,iCell) * dSdzTop(k,iCell)) )
              end do
          end do
          
          ! To do ----
          ! Local estimates of horizontal physical mixing due to laplacian diffusion 
          ! Local estimates of horizontal physical mixing due to biharmonic diffusion 
          ! Local estimates of physical mixing due to redi diffusion 

          ! Local estimates of numerical mixing following Burchard & Rennau (2008) O.M.
          ! Mnum = [ A{c**2} - (A{c})**2 ] / dt_baroclinic 
          ! Order of the mnum tracer tendencies is temp, salt, temp**2,salt**2.

          chiSpurTracerBR08(1,:,:) = ( numericalMixingTracers(3,:,:) - ( numericalMixingTracers(1,:,:)**2 ) )/dt
          chiSpurTracerBR08(2,:,:) = ( numericalMixingTracers(4,:,:) - ( numericalMixingTracers(2,:,:)**2 ) )/dt

          deallocate(dTdzTop)
          deallocate(dSdzTop)
          deallocate(dzTop)

          block => block%next
       end do
 
    end subroutine ocn_compute_discrete_variance_decay!}}}
 
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