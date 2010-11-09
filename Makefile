CC = g++
PROG = bin/gpm_getGNO

# location of GoTools libraries may need to be changed depending on the local computer
GOTLIB	= -L/usr/local/lib/GoTools \
          -lGoTrivariate \
		  -lGoToolsCore

CFLAGS = -Wall -g

# SplineGUI & openGL which should eventually be optionally
GLLIB  = -lglut -lGLU
GUISRC = ../SplineGUI/src/Camera.cpp \
       ../SplineGUI/src/OrthoProjection.cpp \
       ../SplineGUI/src/PointDisplay.cpp \
       ../SplineGUI/src/CurveDisplay.cpp \
       ../SplineGUI/src/SurfaceDisplay.cpp \
       ../SplineGUI/src/VolumeDisplay.cpp \
       ../SplineGUI/src/DisplayObjectSet.cpp \
       ../SplineGUI/src/Button.cpp \
       ../SplineGUI/src/SplineGUI.cpp \
       ../SplineGUI/src/CurvePoint.cpp
GUIINC = ../SplineGUI/include

SRCS   = src/SplineModel.cpp \
         src/TopologySet.cpp \
         src/Line.cpp \
         src/Volume.cpp \
         src/Vertex.cpp \
         src/Face.cpp


LIBS = $(GLLIB) $(GOTLIB) -Iinclude -I$(GUIINC)

all: getGNO gui getPROP refine guiRefine

getGNO: $(SRCS) src/main_getGNO.cpp
	$(CC) $(CFLAGS) -o bin/gpm_getGNO src/main_getGNO.cpp $(SRCS) $(LIBS)

gui: $(SRCS) $(GUISRC) src/main_gui.cpp 
	$(CC) $(CFLAGS) -o bin/gpm_gui src/main_gui.cpp $(SRCS) $(GUISRC) $(LIBS)

getPROP: $(SRCS) src/main_getPROP.cpp 
	$(CC) $(CFLAGS) -o bin/gpm_getPROP src/main_getPROP.cpp $(SRCS) $(LIBS)

refine: $(SRCS) src/main_refine.cpp 
	$(CC) $(CFLAGS) -o bin/gpm_refine src/main_refine.cpp $(SRCS) $(LIBS)

guiRefine: $(SRCS) $(GUISRC) src/main_refine_gui.cpp 
	$(CC) $(CFLAGS) -o bin/gpm_refine_gui src/main_refine_gui.cpp $(SRCS) $(GUISRC) $(LIBS)

clean:
	rm -f bin/gpm_getGNO bin/gpm_gui bin/gpm_getPROP bin/gpm_refine bin/gpm_refine_gui
