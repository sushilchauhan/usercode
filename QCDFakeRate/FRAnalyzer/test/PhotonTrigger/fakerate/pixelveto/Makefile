PACKAGE       = myPlot
ROOTINCLUDES  = -I$(ROOTSYS)/include
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

LIBS          = $(ROOTLIBS) 
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lX11

CC            = gcc

CCFLAGS       = -fPIC -Wall -Wno-deprecated -m64
LDFLAGS       = -fPIC -Wl -m64

LD_LIBRARY_PATH:= $(ROOTSYS)/lib/:$(LD_LIBRARY_PATH)

INCLUDES    = $(ROOTINCLUDES) -I./ 

all:	lib$(PACKAGE).so

$(PACKAGE)dict : $(PACKAGE).h
	@echo "Making $@.C $@.h"
	rm -f $@.C  $@.h
	rootcint $@.C -c $(PACKAGE).h

lib$(PACKAGE).so : $(PACKAGE).C $(PACKAGE)dict
	@echo "Compiling $<"
	$(CC) -shared -o $@ $(CCFLAGS) $(INCLUDES) $(PACKAGE)dict.C $(PACKAGE).C $(LIBS) $(GLIBS)


clean:        
	@echo removing $(PACKAGE)dict.C $(PACKAGE)dict.h $(PACKAGE).so
	rm -f $(PACKAGE)dict.C $(PACKAGE)dict.h lib$(PACKAGE).so
