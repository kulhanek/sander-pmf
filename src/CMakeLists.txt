# ==============================================================================
# SANDER CMake File
# ==============================================================================

INCLUDE_DIRECTORIES(lib/include SYSTEM)
INCLUDE_DIRECTORIES(lib/utils SYSTEM)
INCLUDE_DIRECTORIES(lib/sqm SYSTEM)
INCLUDE_DIRECTORIES(lib/pbsa SYSTEM)
INCLUDE_DIRECTORIES(lib/rism SYSTEM)

ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(bin)
