# ==============================================================================
# SANDER CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(UTILS_SRC
        nxtsec.F
        random.F90
        veclib.F90
        wallclock.c
        AmberNetcdf.F90
        ../pbsa/sa_driver.F90
        ../pbsa/timer.F90
#        mexit.F90
#        amopen.F90
#        assert.F90
#        amflsh.F90
#        sander_bomb.F90
        )

# ------------------------------------------------------------------------------

ADD_DEFINITIONS(${SANDER_DEFS})

# create static library --------------------------------------------------------
IF(LIBS_STATIC)
    ADD_LIBRARY(utils_static STATIC ${UTILS_SRC})
    SET_TARGET_PROPERTIES(utils_static PROPERTIES
                            OUTPUT_NAME utils
                            CLEAN_DIRECT_OUTPUT 1)
ENDIF(LIBS_STATIC)

# create shared library --------------------------------------------------------
IF(LIBS_SHARED)
    ADD_LIBRARY(utils_shared SHARED ${UTILS_SRC})
    SET_TARGET_PROPERTIES(utils_shared PROPERTIES
                            OUTPUT_NAME utils
                            CLEAN_DIRECT_OUTPUT 1
                            VERSION ${SANDER_PMF_VERSION})
                            
	TARGET_LINK_LIBRARIES(utils_shared
	        ${FPMF_LIBRARIES}
	        )                            
                            
    INSTALL(TARGETS utils_shared
            DESTINATION lib)
ENDIF(LIBS_SHARED)

# ------------------------------------------------------------------------------
