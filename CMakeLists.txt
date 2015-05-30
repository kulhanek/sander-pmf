# ==============================================================================
# SANDER 11.0 CMake File - Taken from AMBER 11 distribution archive
# ==============================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
#SET(CMAKE_VERBOSE_MAKEFILE ON)

# ==============================================================================
# conditional setup ------------------------------------------------------------
# ==============================================================================

SET(LIBS_SHARED             ON  CACHE BOOL "Should the shared version of pmflib library be built?")
SET(LIBS_STATIC             OFF CACHE BOOL "Should the static version of pmflib library be built?")

# ==============================================================================
# project setup ----------------------------------------------------------------
# ==============================================================================

# define project name ----------------------------
SET(SANDER_PMF_VERSION 14.04)

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

# BLAS/LAPACK ===============
IF(NOT MKL_FOUND)
    # BLAS ======================
    SET(BLAS_ROOT ${DEVELOPMENT_ROOT}/thirdparty/blas/netlib)
    LINK_DIRECTORIES(${BLAS_ROOT}/lib)
    SET(BLAS_LIB_NAME fblasd)

    # LAPACK =====================
    SET(LAPACK_ROOT ${DEVELOPMENT_ROOT}/thirdparty/lapack/3.1.1)
    LINK_DIRECTORIES(${LAPACK_ROOT}/lib)
    SET(LAPACK_LIB_NAME flapackd)
ELSE(NOT MKL_FOUND)
    SET(LAPACK_LIB_NAME ${MKL_LIBS})
ENDIF(NOT MKL_FOUND)

# NetCDFCore ==================
SET(NETCDF_ROOT ${DEVELOPMENT_ROOT}/thirdparty/netcdfcore/4.0.1)
LINK_DIRECTORIES(${NETCDF_ROOT}/lib)
INCLUDE_DIRECTORIES(${NETCDF_ROOT}/src/lib/fnetcdf SYSTEM)
INCLUDE_DIRECTORIES(${NETCDF_ROOT}/include SYSTEM)
SET(NETCDF_FLIB_NAME fnetcdf)

# FFTW3 =======================
SET(FFTW3_ROOT ${DEVELOPMENT_ROOT}/thirdparty/fftw/3.3.4)
#IF(EXISTS "${FFTW3_ROOT}/")
#    LINK_DIRECTORIES(${FFTW3_ROOT}/lib)
#    SET(FFTW3_LIB_NAME fftw3)
#    SET(HAVE_FFTW3 TRUE)
#ENDIF(EXISTS "${FFTW3_ROOT}/")

# PMFLIB ======================
SET(PMFLIB_ROOT ${DEVELOPMENT_ROOT}/projects/pmflib/4.0)
INCLUDE_DIRECTORIES(${PMFLIB_ROOT}/src/lib/fdrivers/sander SYSTEM)
LINK_DIRECTORIES(${PMFLIB_ROOT}/lib)
LINK_DIRECTORIES(${PMFLIB_ROOT}/lib/drivers)
#SET(PMFLIB_FLIB_NAME fpmf)
#SET(PMFLIB_DRIVER_LIB_NAME fpmfdrv_sander)

# MPI ======================
GET_FILENAME_COMPONENT(FCOMP_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)
if(FCOMP_NAME STREQUAL mpif90)
    ADD_DEFINITIONS(-DMPI)
    SET(MPI_ENABLED TRUE)
    SET(MPI_HOME "$ENV{MPI_HOME}")
    MESSAGE("-- MPI status: enabled")
    MESSAGE("-- MPI home:   ${MPI_HOME}")
ENDIF(FCOMP_NAME STREQUAL mpif90)

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
        sqm pbsa sebomd rism
        )
    SET(SANDER_DEFS -DNDEBUG -DBINTRAJ -DPMFLIB -DRISMSANDER)
ELSE(HAVE_FFTW3)
    SET(FPMF_LIBRARIES
        ${PMFLIB_DRIVER_LIB_NAME} ${PMFLIB_FLIB_NAME} ${NETCDF_FLIB_NAME}
        ${LAPACK_LIB_NAME} ${BLAS_LIB_NAME}
        )

    SET(SANDER_LIBRARIES
        sqm pbsa sebomd
        )
    SET(SANDER_DEFS -DNDEBUG -DBINTRAJ -DPMFLIB)
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

    IF(MPI_HOME)
        SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH}:${MPI_HOME}/lib)
    ENDIF(MPI_HOME)

    SET(INTEL_CORE_HOME $ENV{INTEL_CORE_HOME})
    IF(INTEL_CORE_HOME)
        SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH}:${INTEL_CORE_HOME})
    ENDIF(INTEL_CORE_HOME)

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

    # Fortran options
    SET(CMAKE_Fortran_FLAGS_RELEASE "-O2 -ffree-line-length-none")
    SET(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ffree-line-length-none")

    # fix the problem with dynamic libraries and xlf compiler on SP5
    IF(FCOMP_NAME STREQUAL xlf)
        SET(CMAKE_Fortran_FLAGS_RELEASE "-O2 -qsuppress=1501-510:cmpmsg")
        SET(CMAKE_Fortran_FLAGS_DEBUG   "-g")
        IF(LIBS_SHARED)
            MESSAGE(FATAL_ERROR "LIBS_SHARED is ON, which is not supported with xlf compiler!")
        ENDIF(LIBS_SHARED)
    ELSEIF(FCOMP_NAME STREQUAL ifort)
        # ifort
        SET(CMAKE_Fortran_FLAGS_RELEASE "-O2 -diag-disable 8291 -diag-disable 8290")
        SET(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback -check all -diag-disable 8291 -diag-disable 8290")
    ENDIF(FCOMP_NAME STREQUAL xlf)

    SET(CMAKE_BUILD_TYPE Release)
    #SET(CMAKE_BUILD_TYPE Debug)
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

# ==============================================================================
# project subdirectories  ------------------------------------------------------
# ==============================================================================

ADD_SUBDIRECTORY(src)

IF(MPI_ENABLED)
    INSTALL(FILES
                share/bin/sander.PMF
            DESTINATION
                bin
            PERMISSIONS
                OWNER_EXECUTE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
            )
ENDIF(MPI_ENABLED)

