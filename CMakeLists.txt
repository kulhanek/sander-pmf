# ==============================================================================
# SANDER 17.0 CMake File - Taken from AMBER16/AmberTools17 distribution archive
# ==============================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 3.10 FATAL_ERROR)
#SET(CMAKE_VERBOSE_MAKEFILE ON)
PROJECT(SANDER-PMF)

# ==============================================================================
# conditional setup ------------------------------------------------------------
# ==============================================================================

SET(LIBS_SHARED             ON  CACHE BOOL "Should the shared version of pmflib library be built?")
SET(LIBS_STATIC             OFF CACHE BOOL "Should the static version of pmflib library be built?")

# ==============================================================================
# project setup ----------------------------------------------------------------
# ==============================================================================

# define project name ----------------------------
SET(SANDER_PMF_VERSION 17.06)

# what languages we use --------------------------
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(Fortran)

# ==============================================================================
# external dependencies --------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    SET(DEVELOPMENT_ROOT "$ENV{DEVELOPMENT_ROOT}")
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

IF(NOT DEVELOPMENT_ROOT)
    MESSAGE(FATAL_ERROR "DEVELOPMENT_ROOT variable is not set!")
ENDIF(NOT DEVELOPMENT_ROOT)

# ==============================================================================
# external libraries -----------------------------------------------------------
# ==============================================================================
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${DEVELOPMENT_ROOT}/projects/cmake-modules/1.0/share/modules")

# setup IntelCDK if available
INCLUDE(IntelCDK)

# setup MPI if available
INCLUDE(OpenMPI)

# boost =======================
FIND_PACKAGE(MyBoost REQUIRED)
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})

# BLAS/LAPACK =================
# setup BLAS/LAPACK either from MKL or from fallback netlib blas and lapack
FIND_PACKAGE(BlasLapack)

# NetCDFCore ==================
SET(NETCDF_ROOT ${DEVELOPMENT_ROOT}/thirdparty/netcdfcore/4.0.1)
LINK_DIRECTORIES(${NETCDF_ROOT}/lib)
INCLUDE_DIRECTORIES(${NETCDF_ROOT}/src/lib/fnetcdf SYSTEM)
INCLUDE_DIRECTORIES(${NETCDF_ROOT}/include SYSTEM)
SET(NETCDF_FLIB_NAME fnetcdf)

## FFTW3 =======================
#SET(FFTW3_ROOT ${DEVELOPMENT_ROOT}/thirdparty/fftw/3.3.4)
#IF(EXISTS "${FFTW3_ROOT}/")
#    LINK_DIRECTORIES(${FFTW3_ROOT}/lib)
#    INCLUDE_DIRECTORIES(${FFTW3_ROOT}/include SYSTEM)
#    SET(FFTW3_LIB_NAME fftw3)
#    IF(MPI_ENABLED)
#        SET(FFTW3_LIB_NAME ${FFTW3_LIB_NAME} fftw3_mpi)     
#    ENDIF(MPI_ENABLED)
#    SET(HAVE_FFTW3 TRUE)
#ENDIF(EXISTS "${FFTW3_ROOT}/")

# PMFLIB ======================
SET(PMFLIB_ROOT ${DEVELOPMENT_ROOT}/projects/pmflib/6.0)
INCLUDE_DIRECTORIES(${PMFLIB_ROOT}/src/lib/fdrivers/sander SYSTEM)
LINK_DIRECTORIES(${PMFLIB_ROOT}/lib)
LINK_DIRECTORIES(${PMFLIB_ROOT}/lib/drivers)
SET(PMFLIB_FLIB_NAME fpmf)
SET(PMFLIB_DRIVER_LIB_NAME fpmfdrv_sander)

# policies ---------------------------------------
IF(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

# ==============================================================================

IF(HAVE_FFTW3)
    SET(FPMF_LIBRARIES
        ${PMFLIB_DRIVER_LIB_NAME} ${PMFLIB_FLIB_NAME} ${NETCDF_FLIB_NAME}
        ${LAPACK_LIB_NAME} ${BLAS_LIB_NAME} ${FFTW3_LIB_NAME}
        )

    SET(SANDER_LIBRARIES
        sqm pbsa sebomd rism utils
        )
    SET(SANDER_DEFS -DNDEBUG -DBINTRAJ -DPMFLIB -DRISMSANDER -DDISABLE_NFE) 
ELSE(HAVE_FFTW3)
    SET(FPMF_LIBRARIES
        ${PMFLIB_DRIVER_LIB_NAME} ${PMFLIB_FLIB_NAME} ${NETCDF_FLIB_NAME}
        ${LAPACK_LIB_NAME} ${BLAS_LIB_NAME}
        )

    SET(SANDER_LIBRARIES
        sqm pbsa sebomd utils
        )
    SET(SANDER_DEFS -DNDEBUG -DBINTRAJ -DPMFLIB -DDISABLE_NFE)
ENDIF(HAVE_FFTW3)

# ==============================================================================
# output directories -----------------------------------------------------------
# ==============================================================================

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)

# ==============================================================================
# rpaths -----------------------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

    # the RPATH to be used when installing
    SET(CMAKE_INSTALL_RPATH "${FFTW3_ROOT}/lib:${BLAS_ROOT}/lib:${LAPACK_ROOT}/lib:${NETCDF_ROOT}/lib:${PMFLIB_ROOT}/lib:${PMFLIB_ROOT}/lib/drivers:${CMAKE_BINARY_DIR}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

# ==============================================================================
# compiler setups --------------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    # C options
    SET(CMAKE_C_FLAGS_RELEASE "-O2")
    SET(CMAKE_C_FLAGS_DEBUG "-O0 -g")

    # get the name of fortran compiler
    get_filename_component(FCOMP_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)
    IF(MPI_ENABLED)
        EXECUTE_PROCESS(COMMAND mpifort --showme:command OUTPUT_VARIABLE FCOMP_NAME OUTPUT_STRIP_TRAILING_WHITESPACE)
    ENDIF(MPI_ENABLED)       

    # Fortran options
    SET(CMAKE_Fortran_FLAGS_RELEASE "-O2 -ffree-line-length-none")
    SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ffree-line-length-none")

    # fix the problem with dynamic libraries and xlf compiler on SP5
    IF(FCOMP_NAME STREQUAL ifort)
        # ifort
        SET(CMAKE_Fortran_FLAGS_RELEASE "-O2 -diag-disable 8291 -diag-disable 8290")
        SET(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback -check all -diag-disable 8291 -diag-disable 8290")
    ENDIF(FCOMP_NAME STREQUAL ifort)

    SET(CMAKE_BUILD_TYPE Release)
    #SET(CMAKE_BUILD_TYPE Debug)
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

# ==============================================================================
# project subdirectories  ------------------------------------------------------
# ==============================================================================

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(share)

