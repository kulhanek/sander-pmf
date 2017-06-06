# ==============================================================================
# SANDER CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(PBSA_SRC
#        sa_driver.F90
        density.F90
        aug_iccg.F90
        bicg.F90
        coed20.F90
        coed6.F90
        curv.F90
        d1mach.F90
        daxpy.F90
        dbcg.F90
        dchkw.F90
        dcopy.F90
        ddot.F90
        debug.F90
        decomp.F90
        dgmres.F90
        dhels.F90
        dheqr.F90
        dnrm2.F90
        dorth.F90
        dpigmr.F90
        drlcal.F90
        ds2y.F90
        dscal.F90
        dsilus.F90
        dslubc.F90
        dslugm.F90
        dslui2.F90
        dslui4.F90
        dslui.F90
        dsluti.F90
        dsmtv.F90
        dsmv.F90
        dsvdc.F90
        dxlcal.F90
        egb.F90
        ene.F90
        fdump.F90
        fftw3.F90
        force.F90
        gen_dx_file.F90
        getcoor.F90
        gmres.F90
        gmresX.F90
        GrToPr.F90
        i1mach.F90
        IIM.F90
        iimod.F90
        indexg.F90
        interpX.F90
        irre31.F90
        irre32.F90
        isdbcg.F90
        isdgmr.F90
        j4save.F90
        jumps.F90
        locmem.F90
        matvec3.F90
#        membrane.F90
        memory_module.F90
        miniop.F90
        myopen.F90
        np_force.F90
        parms.F90
        pb_augdrv.F90
        pb_chkfrc.F90
        pb_direct.F90
#        pb_exmol.F90
        pb_fddrv.F90
#        pb_fft.F90
#        pb_fftsolv.F90
        pb_force.F90
        pb_iimdrv.F90
        pb_init.F90
        pb_list.F90
        pb_lsolver.F90
        pb_mpfrc.F90
        pb_nlsolver.F90
        pb_read.F90
        pbsaapi.F90
#        pbsa.F90
        pbsa_lib.F90
        pb_write.F90
        phidev.F90
        problem.F90
        prodis.F90
        project.F90
        pythag.F90
        qint.F90
        qld.F90
        qs2i1d.F90
        rdparm.F90
        rdpqr.F90
        regular.F90
        rfree.F90
        rgroup.F90
        runmd.F90
        runmin.F90
        svbksb.F90
        svdcmp.F90
        timer.F90
        transf.F90
        wint.F90
        xercnt.F90
        xerhlt.F90
        xermsg.F90
        xerprn.F90
        xersve.F90
        xgetua.F90
        pb_fdfrc.F90
        pb_p3m.F90
        pb_bldsys.F90
        pb_crgview.F90
#        random.F90
        )

# ------------------------------------------------------------------------------

ADD_DEFINITIONS(${SANDER_DEFS})

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(pbsa_static STATIC ${PBSA_SRC})
    SET_TARGET_PROPERTIES(pbsa_static PROPERTIES
                            OUTPUT_NAME pbsa
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(pbsa_shared SHARED ${PBSA_SRC})
    ADD_DEPENDENCIES(pbsa_shared utils_shared)

    SET_TARGET_PROPERTIES(pbsa_shared PROPERTIES
                            OUTPUT_NAME pbsa
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${SANDER_PMF_VERSION})
	TARGET_LINK_LIBRARIES(pbsa_shared
			utils
	        ${FPMF_LIBRARIES}
	        )                            
    INSTALL(TARGETS pbsa_shared
            DESTINATION lib)
ENDIF(LIBS_SHARED)

# ------------------------------------------------------------------------------
