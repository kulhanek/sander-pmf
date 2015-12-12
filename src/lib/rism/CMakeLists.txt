# ==============================================================================
# SANDER CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(RISM_SRC
        amber_rism_interface.F90
        array_util.F90
        blend.F90
        constants.F90
        erfcfun.F90
        fce_c.F90
        fftw3.F90
        getopts_c.F90
        mdiis_blas2_c.F90
        mdiis_blas_c.F90
        mdiis_c.F90
        mdiis_orig_c.F90
        pubfft.F90
        quaternion.F90
        rism1d_c.F90
        rism1d_closure_c.F90
#        rism1d.F90                                                                                                                                                    
        rism1d_hnc_c.F90                                                                                                                                              
        rism1d_kh_c.F90
        rism1d_mv0_c.F90
        rism1d_potential_c.F90
        rism1d_psen_c.F90
        rism1d_py_c.F90
        rism3d_c.F90
        rism3d_closure_c.F90
        rism3d_csv.F90
        rism3d_debug_c.F90
        rism3d_fft.F90
        rism3d_grid_c.F90
        rism3d_hnc_c.F90
        rism3d_kh_c.F90
        rism3d_opendx.F90
        rism3d.orave.F90
        rism3d_potential_c.F90
        rism3d_psen_c.F90
        rism3d_solu_c.F90
        rism3d_solv_c.F90
#        rism3d.thermo.F90
        rism3d_xyzv.F90
        rism_parm.F90
        rism_report_c.F90
        rism_timer_c.F90
        rism_util.F90
        safemem.F90
        solvmdl_c.F90
#        volslice.F90

        xblas/BLAS_daxpby_x-f2c.c
        xblas/BLAS_daxpby_x.c
        xblas/BLAS_dwaxpby_x-f2c.c
        xblas/BLAS_dwaxpby_x.c
        xblas/BLAS_dgemm_x-f2c.c
        xblas/BLAS_dgemm_x.c
        xblas/BLAS_error.c
        )

# ------------------------------------------------------------------------------

ADD_DEFINITIONS(${SANDER_DEFS} -DSANDER -DCONFIG_FC_UNDERSCORE)

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(rism_static STATIC ${RISM_SRC})
    SET_TARGET_PROPERTIES(sqm_static PROPERTIES
                            OUTPUT_NAME rism
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(rism_shared SHARED ${RISM_SRC})
    ADD_DEPENDENCIES(rism_shared utils_shared)
    SET_TARGET_PROPERTIES(rism_shared PROPERTIES
                            OUTPUT_NAME rism
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${SANDER_PMF_VERSION})
    INSTALL(TARGETS rism_shared
            DESTINATION lib)

ENDIF(LIBS_SHARED)

# ------------------------------------------------------------------------------
