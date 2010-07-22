CC = g++
PROG = bin/gpm

# location of GoTools libraries may need to be changed depending on the local computer
GOTLIB	= -L/usr/local/lib/GoTools \
          -lGoTrivariate \
		  -lGoToolsCore

CFLAGS = -Wall -g

SRCS =	src/main.cpp \
      	src/SplineModel.cpp \
      	src/TopologySet.cpp \
      	src/Line.cpp \
      	src/Volume.cpp \
      	src/Face.cpp


LIBS = $(GOTLIB) -Iinclude

all: $(PROG)

$(PROG):	$(SRCS)
	$(CC) $(CFLAGS) -o $(PROG) $(SRCS) $(LIBS)

clean:
	rm -f $(PROG)
