# ==============================================================================
# sander-pmf
# ==============================================================================

IF(MPI_ENABLED)
    CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/bin/sander.PMF.in ${CMAKE_CURRENT_SOURCE_DIR}/bin/sander.PMF @ONLY)
    INSTALL(FILES
                bin/sander.PMF
            DESTINATION
                bin
            PERMISSIONS
                OWNER_EXECUTE OWNER_READ OWNER_WRITE GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
            )
ENDIF(MPI_ENABLED)
