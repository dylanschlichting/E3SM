include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

#message(STATUS "gcp PROJECT_NAME=${PROJECT_NAME} USE_CUDA=${USE_CUDA} KOKKOS_ENABLE_CUDA=${KOKKOS_ENABLE_CUDA}")
# use default backend?

if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if (BUILD_THREADED)
    include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
  else()
    include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
  endif()
else()
  include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
endif()

include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)

set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

#message(STATUS "gcp CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()
