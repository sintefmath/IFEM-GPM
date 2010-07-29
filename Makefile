CC = g++
PROG = bin/gpm

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

SRCS =	src/main_fenris.cpp \
      	src/SplineModel.cpp \
      	src/TopologySet.cpp \
      	src/Line.cpp \
      	src/Volume.cpp \
      	src/Vertex.cpp \
      	src/Face.cpp \
		$(FENRISSRC)


LIBS = $(GLLIB) $(GOTLIB) -Iinclude -I$(FENRISINC)

all: $(PROG)

$(PROG):	$(SRCS)
	$(CC) $(CFLAGS) -o $(PROG) $(SRCS) $(LIBS)

clean:
	rm -f $(PROG)
