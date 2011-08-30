PROJECT(GPM)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ENABLE_LANGUAGE(CXX)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${CMAKE_SOURCE_DIR}/cmake/Modules)

# Required packages
FIND_PACKAGE(GoTools REQUIRED)
FIND_PACKAGE(GoTrivariate REQUIRED)
FIND_PACKAGE(GLUT REQUIRED)

# Required libraries
SET(DEPLIBS ${GoTrivariate_LIBRARIES} ${GoTools_LIBRARIES}
            ${GLUT_glut_LIBRARY}
)

# Required include directories
SET(INCLUDES
  ${GoTools_INCLUDE_DIRS}
  ${GoTrivariate_INCLUDE_DIRS}
  ${GLUT_INCLUDE_DIR}
  ${PROJECT_SOURCE_DIR}/include
)

INCLUDE_DIRECTORIES(${INCLUDES})

SET(EXECUTABLE_OUTPUT_PATH bin)
SET(LIBRARY_OUTPUT_PATH lib)

# This line may be needed only on AMD 64 bit... @jbt
# 091120: Adding this because GUV-plugins will need it. (J.O.)
IF(NOT WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF(NOT WIN32)

# Enable all warnings
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# Make the GPM library
FILE(GLOB_RECURSE GPM_SRCS ${PROJECT_SOURCE_DIR}/src/SplineModel.cpp
                           ${PROJECT_SOURCE_DIR}/src/TopologySet.cpp
                           ${PROJECT_SOURCE_DIR}/src/Line.cpp
                           ${PROJECT_SOURCE_DIR}/src/Volume.cpp
                           ${PROJECT_SOURCE_DIR}/src/Vertex.cpp
                           ${PROJECT_SOURCE_DIR}/src/Face.cpp)
ADD_LIBRARY(GPM ${GPM_SRCS})
TARGET_LINK_LIBRARIES(GPM ${DEPLIBS})

# Make the Apps
FILE(GLOB_RECURSE getGNO_SRCS ${PROJECT_SOURCE_DIR}/src/main_getGNO.cpp)
ADD_EXECUTABLE(getGNO ${getGNO_SRCS})
TARGET_LINK_LIBRARIES(getGNO GPM ${DEPLIBS})

FILE(GLOB_RECURSE getPROP_SRCS ${PROJECT_SOURCE_DIR}/src/main_getPROP.cpp)
ADD_EXECUTABLE(getPROP ${getPROP_SRCS})
TARGET_LINK_LIBRARIES(getPROP GPM ${DEPLIBS})

FILE(GLOB_RECURSE refine_SRCS ${PROJECT_SOURCE_DIR}/src/main_refine.cpp)
ADD_EXECUTABLE(refine ${refine_SRCS})
TARGET_LINK_LIBRARIES(refine GPM ${DEPLIBS})

FILE(GLOB_RECURSE orderElevation_SRCS ${PROJECT_SOURCE_DIR}/src/main_orderElevation.cpp)
ADD_EXECUTABLE(orderElevation ${orderElevation_SRCS})
TARGET_LINK_LIBRARIES(orderElevation ${DEPLIBS})

IF(NOT GPM_INSTALL_PREFIX)
  SET(GPM_INSTALL_PREFIX $ENV{HOME})
ENDIF(NOT GPM_INSTALL_PREFIX)

# 'install' target
IF(WIN32)
  # TODO
ELSE(WIN32)
  # lib
  INSTALL(TARGETS GPM DESTINATION ${GPM_INSTALL_PREFIX}/lib)

  # headers
  FILE(GLOB_RECURSE GPM_HEADERS include/*.h)
  INSTALL(FILES ${GPM_HEADERS}
          DESTINATION ${GPM_INSTALL_PREFIX}/include/GPM)

  FILE(GLOB_RECURSE LOCAL_CMAKE_MODULES cmake/Modules/*.cmake)

  # cmake modules
  INSTALL(FILES ${LOCAL_CMAKE_MODULES} DESTINATION ${GPM_INSTALL_PREFIX}/cmake/Modules)
ENDIF(WIN32)
