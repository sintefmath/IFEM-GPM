# - Tries to find the GoTools Core library
#
# Written by: jan.b.thomassen@sintef.no
#


# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(GoTools_INCLUDE_DIRS ${GoToolsCore_SOURCE_DIR}/include
    CACHE PATH "Path to GoTools Core header files")
  # Library
  SET(GoTools_LIBRARIES GoToolsCore
    CACHE FILE "GoTools Core library")
ENDIF(GoTools_BUILD_ALL)

# Find header files
FIND_PATH(GoTools_INCLUDE_DIRS "GoTools/geometry/SplineSurface.h"
  /sima/libs/GoTools/include
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  "C:/Program Files (x86)/GoTools/include"
  )


# Find library
FIND_LIBRARY(GoTools_LIBRARIES
  NAMES GoToolsCore 
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "C:/Program Files (x86)/GoTools/lib"
  /sima/libs/GoTools/lib
  PATH_SUFFIXES GoTools
  )

# Check for newer GoTools
EXECUTE_PROCESS(COMMAND cat "${GoTools_INCLUDE_DIRS}/GoTools/geometry/GoTools.h" OUTPUT_VARIABLE GOTOOLS_HEADER)
STRING(REGEX REPLACE ".*GO_VERSION_MAJOR ([0-9]+).*" "\\1" GoTools_VERSION_MAJOR ${GOTOOLS_HEADER})
STRING(REGEX REPLACE ".*GO_VERSION_MINOR ([0-9]+).*" "\\1" GoTools_VERSION_MINOR ${GOTOOLS_HEADER})

IF ("${GoTools_VERSION_MAJOR}" MATCHES 3)
  INCLUDE(CheckCXXCompilerFlag)
  IF(CMAKE_CXX_COMPILER_ID MATCHES GNU)
  # check if compiler supports c++-0x
    CHECK_CXX_COMPILER_FLAG("-std=gnu++0x" HAVE_0x)
    IF(HAVE_0x)
      SET(GoTools_CXX_FLAGS "-std=gnu++0x")
    ELSE(HAVE_0x)
      MESSAGE("A compiler with c++-0x support is needed")
      EXIT(1)
    ENDIF(HAVE_0x)
  ELSE(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    MESSAGE("Cannot verify that compiler supports c++-0x, bailing")
    EXIT(1)
  ENDIF(CMAKE_CXX_COMPILER_ID MATCHES GNU)
ENDIF ("${GoTools_VERSION_MAJOR}" MATCHES 3)

# Check that we have found everything
SET(GoTools_FOUND FALSE)
IF(GoTools_INCLUDE_DIRS AND GoTools_LIBRARIES)
  SET(GoTools_FOUND TRUE)
ENDIF(GoTools_INCLUDE_DIRS AND GoTools_LIBRARIES)
