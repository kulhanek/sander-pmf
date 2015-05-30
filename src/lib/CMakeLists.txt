# ==============================================================================
# SANDER CMake File
# ==============================================================================

ADD_SUBDIRECTORY(sqm)
ADD_SUBDIRECTORY(sebomd)
ADD_SUBDIRECTORY(pbsa)
IF(HAVE_FFTW3)
    ADD_SUBDIRECTORY(rism)
ENDIF(HAVE_FFTW3)
