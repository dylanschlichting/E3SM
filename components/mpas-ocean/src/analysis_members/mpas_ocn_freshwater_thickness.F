! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_freshwater_thickness
!
!> \brief MPAS ocean analysis mode member: Computes freshwater quantities
!>        from the riverTagTracers.

!> \author Dylan Schlichting
!> \date   February 06, 2025, updated March 26, 2025
!> \details
!>   Computes freshwater thickness by depth-integrating river tracers
!>   Details described in Zhang et al. 2012 JGR Oceans
!>   https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JC008108
!>   Returns a nCells x Time array for each river.
!> 
!>   fw_thickness = \int river_dye * dz
!>   fw_volume = \iiint river_dye * dz
!>   fw_u_flux = \int river_dye * u * dz
!>   fw_v_flux = \int river_dye * v * dz
!>   The fluxes are used for plotting vectors on plan view plots of thickness.
!-----------------------------------------------------------------------

module ocn_freshwater_thickness

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_variables

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

   public :: ocn_init_freshwater_thickness, &
             ocn_compute_freshwater_thickness, &
             ocn_restart_freshwater_thickness, &
             ocn_finalize_freshwater_thickness

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_freshwater_thickness
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author Dylan Schlichting
!> \date    February 06, 2025, updated March 26, 2025
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_freshwater_thickness(domain, err)!{{{

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

   end subroutine ocn_init_freshwater_thickness!}}}

!***********************************************************************
!
!  routine ocn_compute_freshwater_thickness
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Dylan Schlichting
!> \date    February 06, 2025
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_freshwater_thickness(domain, timeLevel, err)!{{{

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

      type(mpas_pool_type), pointer :: freshwaterThicknessAMPool
      type(dm_info) :: dminfo
      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: statePool
      type(mpas_pool_type), pointer :: meshPool
      type(mpas_pool_type), pointer :: tracersPool

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve
      integer :: iTracer, k, iCell
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot

      real(kind=RKIND), dimension(:, :, :), pointer :: tracers, riverTagTracers
      real(kind=RKIND), dimension(:, :), pointer :: riverTagTracer1RiverRunoff, riverTagTracer2RiverRunoff
      real(kind=RKIND), dimension(:, :), pointer :: layerThickness, velocityZonal, velocityMeridional
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real(kind=RKIND), dimension(:), pointer ::  riverTagTracer1FwThickness, riverTagTracer2FwThickness, & 
                                                  riverTagTracer1FwVolume, riverTagTracer2FwVolume, &
                                                  riverTagTracer1ZonalFwFlux, riverTagTracer1MeridionalFwFlux, &
                                                  riverTagTracer2ZonalFwFlux, riverTagTracer2MeridionalFwFlux

      real(kind=RKIND) :: freshwaterLayer_river1, freshwaterLayer_river2, &
                          freshwaterLayerVolume_river1, freshwaterLayerVolume_river2, &
                          zonalFreshwaterFluxLayer_river1, meridionalFreshwaterFluxLayer_river1, & 
                          zonalFreshwaterFluxLayer_river2, meridionalFreshwaterFluxLayer_river2

      err = 0

      dminfo = domain%dminfo

      block => domain%blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block%structs, 'state', statePool)
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block%structs, 'freshwaterThickness', freshwaterThicknessAMPool)
         call mpas_pool_get_dimension(block%dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block%dimensions, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(tracersPool, 'riverTagTracers', riverTagTracers, timeLevel)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
         call mpas_pool_get_array(freshwaterThicknessAMPool, 'riverTagTracer1FwThickness', riverTagTracer1FwThickness) 
         call mpas_pool_get_array(freshwaterThicknessAMPool, 'riverTagTracer2FwThickness', riverTagTracer2FwThickness) 
         call mpas_pool_get_array(freshwaterThicknessAMPool, 'riverTagTracer1FwVolume', riverTagTracer1FwVolume) 
         call mpas_pool_get_array(freshwaterThicknessAMPool, 'riverTagTracer2FwVolume', riverTagTracer2FwVolume) 
         call mpas_pool_get_array(freshwaterThicknessAMPool, 'riverTagTracer1ZonalFwFlux', riverTagTracer1ZonalFwFlux) 
         call mpas_pool_get_array(freshwaterThicknessAMPool, 'riverTagTracer1MeridionalFwFlux', riverTagTracer1MeridionalFwFlux) 
         call mpas_pool_get_array(freshwaterThicknessAMPool, 'riverTagTracer2ZonalFwFlux', riverTagTracer2ZonalFwFlux) 
         call mpas_pool_get_array(freshwaterThicknessAMPool, 'riverTagTracer2MeridionalFwFlux', riverTagTracer2MeridionalFwFlux) 

         do iCell = 1, nCellsSolve
            k = 1
            riverTagTracer1FwThickness(iCell) = 0.0_RKIND
            riverTagTracer2FwThickness(iCell) = 0.0_RKIND
            riverTagTracer1FwVolume(iCell) = 0.0_RKIND
            riverTagTracer2FwVolume(iCell) = 0.0_RKIND
            riverTagTracer1ZonalFwFlux(iCell) = 0.0_RKIND
            riverTagTracer1MeridionalFwFlux(iCell) = 0.0_RKIND
            riverTagTracer2ZonalFwFlux(iCell) = 0.0_RKIND
            riverTagTracer2MeridionalFwFlux(iCell) = 0.0_RKIND
            do while (k .le. maxLevelCell(iCell))
               ! Freshwater concentration in each layer > dye*dz
               freshwaterLayer_river1 = layerThickness(k, iCell)*riverTagTracers(1, k, iCell)
               freshwaterLayer_river2 = layerThickness(k, iCell)*riverTagTracers(2, k, iCell)

               ! Freshwater thickness > \int dye*dz   
               riverTagTracer1FwThickness(iCell) = riverTagTracer1FwThickness(iCell) + freshwaterLayer_river1
               riverTagTracer2FwThickness(iCell) = riverTagTracer2FwThickness(iCell) + freshwaterLayer_river2
               
               ! Freshwater volume in each layer > dye*dV, where dV = dA*dz
               freshwaterLayerVolume_river1 = layerThickness(k, iCell)*areaCell(iCell)*riverTagTracers(1, k, iCell)
               freshwaterLayerVolume_river2 = layerThickness(k, iCell)*areaCell(iCell)*riverTagTracers(2, k, iCell)
               
               ! Freshwater volume > \iiint dye*dV
               riverTagTracer1FwVolume(iCell) = riverTagTracer1FwVolume(iCell)+freshwaterLayerVolume_river1
               riverTagTracer2FwVolume(iCell) = riverTagTracer2FwVolume(iCell)+freshwaterLayerVolume_river2
               
               ! Freshwater fluxes in each layer > dye*u*dz, dye*v*dz
               zonalFreshwaterFluxLayer_river1 = layerThickness(k, iCell)*riverTagTracers(1, k, iCell)*velocityZonal(k,iCell)
               meridionalFreshwaterFluxLayer_river1 = layerThickness(k, iCell)*riverTagTracers(1, k, iCell)*velocityMeridional(k,iCell)

               zonalFreshwaterFluxLayer_river2 = layerThickness(k, iCell)*riverTagTracers(2, k, iCell)*velocityZonal(k,iCell)
               meridionalFreshwaterFluxLayer_river2 = layerThickness(k, iCell)*riverTagTracers(2, k, iCell)*velocityMeridional(k,iCell)
               
               ! Freshwater fluxes> \int dye*u*dz, \int dye*v*dz
               riverTagTracer1ZonalFwFlux(iCell) = riverTagTracer1ZonalFwFlux(iCell)+zonalFreshwaterFluxLayer_river1
               riverTagTracer1MeridionalFwFlux(iCell) = riverTagTracer1MeridionalFwFlux(iCell)+meridionalFreshwaterFluxLayer_river1

               riverTagTracer2ZonalFwFlux(iCell) = riverTagTracer2ZonalFwFlux(iCell)+zonalFreshwaterFluxLayer_river2
               riverTagTracer2MeridionalFwFlux(iCell) = riverTagTracer2MeridionalFwFlux(iCell)+meridionalFreshwaterFluxLayer_river2

               k = k + 1
            enddo
         end do

         block => block%next
      end do

   end subroutine ocn_compute_freshwater_thickness!}}}

!***********************************************************************
!
!  routine ocn_restart_freshwater_thickness
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Dylan Schlichting
!> \date    February 06, 2025
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_freshwater_thickness(domain, err)!{{{

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

   end subroutine ocn_restart_freshwater_thickness!}}}

!***********************************************************************
!
!  routine ocn_finalize_freshwater_thickness
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Dylan Schlichting
!> \date    February 06, 2025
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_freshwater_thickness(domain, err)!{{{

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

   end subroutine ocn_finalize_freshwater_thickness!}}}

end module ocn_freshwater_thickness

! vim: foldmethod=marker
