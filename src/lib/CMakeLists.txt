# ==============================================================================
# SANDER CMake File
# ==============================================================================

ADD_SUBDIRECTORY(pbsa)
ADD_SUBDIRECTORY(sqm)
ADD_SUBDIRECTORY(sebomd)
IF(HAVE_FFTW3)
    ADD_SUBDIRECTORY(rism)
ENDIF(HAVE_FFTW3)
