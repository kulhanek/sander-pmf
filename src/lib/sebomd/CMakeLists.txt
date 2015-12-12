# ==============================================================================
# SANDER CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(SEBOMD_SRC
        block.F
        se_atmchg.F
        se_balance.F
#        sebomd.F
        se_bpair.F
        se_bsort.F
        se_bspline.F
        se_clustsub.F
#        se_corrections.F
        se_debug.F
        se_denful.F
        se_densub.F
        se_diag.F
        se_diagp.F
        se_diat.F
        se_dihedrpbc.F
        se_dipole.F
        se_dmxinter.F
        se_doferm.F
#        se_doscf.F
#        se_edriver.F
        se_elecsr.F
#        se_energy.F
        se_escf.F
        se_esqr.F
        se_etimer.F
        se_fock.F
        se_fockpme.F
        se_fourier.F
        se_fshift.F
        se_gcart.F
        se_gdiat.F
        se_gensub.F
        se_getpar.F
        se_glbpnt.F
        se_goverlp.F
        se_grepul.F
        se_ijpair.F
        se_initp.F
        se_mosub.F
        se_mpi_allreduce.F
        se_opnfil.F
        se_overlp.F
        se_pbcgcres.F
        se_pbcxyz.F
        se_pme_calcb.F
        se_pme_calcchg.F
        se_pme_calccm1.F
        se_pme_calcq.F
        se_pme_calctheta.F
        se_pme_derec.F
        se_pme_direct2.F
        se_pme_qm_direct.F
        se_pme_qm_recip.F
        se_pme_recip.F
#        se_pme_setup.F
        se_pmix.F
        se_rdall.F
        se_rdelem.F
        se_rdnum.F
        se_rdtail.F
        se_repul.F
        se_resclustsub.F
        se_setbox.F
        se_setunit.F
        se_setup.F
        se_sproc.F
        se_whatis.F
        se_wrtdmx.F
        )

# ------------------------------------------------------------------------------

ADD_DEFINITIONS(${SANDER_DEFS})

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(sebomd_static STATIC ${SEBOMD_SRC})
    SET_TARGET_PROPERTIES(sqm_static PROPERTIES
                            OUTPUT_NAME sebomd
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(sebomd_shared SHARED ${SEBOMD_SRC})
    SET_TARGET_PROPERTIES(sebomd_shared PROPERTIES
                            OUTPUT_NAME sebomd
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${SANDER_PMF_VERSION})
	TARGET_LINK_LIBRARIES(sebomd_shared
			utils
	        ${FPMF_LIBRARIES}
	        )                             
    INSTALL(TARGETS sebomd_shared
            DESTINATION lib)
ENDIF(LIBS_SHARED)

# ------------------------------------------------------------------------------
