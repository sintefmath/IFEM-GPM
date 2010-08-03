CC = g++
PROG = bin/gpm_getGNO

# location of GoTools libraries may need to be changed depending on the local computer
GOTLIB	= -L/usr/local/lib/GoTools \
          -lGoTrivariate \
		  -lGoToolsCore

CFLAGS = -Wall -g

# Fenris & openGL which should eventually be optionally
GLLIB      = -lglut
FENRISSRC = ../Fenris/src/Camera.cpp \
       ../Fenris/src/OrthoProjection.cpp \
	   ../Fenris/src/PointDisplay.cpp \
	   ../Fenris/src/CurveDisplay.cpp \
	   ../Fenris/src/SurfaceDisplay.cpp \
	   ../Fenris/src/VolumeDisplay.cpp \
	   ../Fenris/src/DisplayObjectSet.cpp \
	   ../Fenris/src/Button.cpp \
	   ../Fenris/src/Fenris.cpp \
	   ../Fenris/src/CurvePoint.cpp
FENRISINC = ../Fenris/include

SRCS =	src/SplineModel.cpp \
      	src/TopologySet.cpp \
      	src/Line.cpp \
      	src/Volume.cpp \
      	src/Vertex.cpp \
      	src/Face.cpp \
		$(FENRISSRC)


LIBS = $(GLLIB) $(GOTLIB) -Iinclude -I$(FENRISINC)

all: getGNO gui getPROP

getGNO: $(SRCS) src/main_getGNO.cpp
	$(CC) $(CFLAGS) -o bin/gpm_getGNO src/main_getGNO.cpp $(SRCS) $(LIBS)

gui: $(SRCS) src/main_gui.cpp 
	$(CC) $(CFLAGS) -o bin/gpm_gui src/main_gui.cpp $(SRCS) $(LIBS)

getPROP: $(SRCS) src/main_getPROP.cpp 
	$(CC) $(CFLAGS) -o bin/gpm_getPROP src/main_getPROP.cpp $(SRCS) $(LIBS)

clean:
	rm -f src/gpm_getGNO src/gpm_gui src/gpm_getPROP
