.IGNORE:

SHELL = /bin/sh

prefix = /home/arman/local
BBH_SYSTEM = LINUX_INTEL64
bindir = $(prefix)/bin

BBH_SYSTEM = LINUX_INTEL64

RNPL   = rnpl

AR     = xiar
RANLIB = ranlib


LIBS       = -lm   -lpamrw  -lamrdw  -lbbhutil  -lvutil -lsvml -lifcore -lm
LDFLAGS    = -O3    -L/global/software/intel/fc/11.0.074/lib/intel64 -L/global/software/openmpi-1.4.4/intel-2011/lib -L/usr/X11R6/lib -L.    -L/home/arman/local/lib

CC       = /global/software/openmpi-1.4.4/intel-2011/bin/mpicc
CFLAGS     = -O3  
CPPFLAGS = -I. -I/global/software/openmpi-1.4.4/intel-2011/include -wd177,266 -DINTEL_8_OR_LATER -DLINUX -DPORT_GROUP -DWant_c_files -I$(prefix)/include   -I/home/arman/local/include

CC_COMP  = $(CC) -c $(CFLAGS) $(CPPFLAGS)
CC_LOAD  = $(CC) $(LDFLAGS)

F77      = /global/software/openmpi-1.4.4/intel-2011/bin/mpif77
F77FLAGS = -O3  -w -w90 -w95  -cm -Vaxlib

F77_COMP   = $(F77) -c $(F77FLAGS)
F77_LOAD   = $(F77) $(F77FLAGS) $(F77_LDFLAGS) $(LDFLAGS)

EXECUTABLES     = wave
all: $(EXECUTABLES)

install: all
	if test ! -d $(prefix); then\
		mkdir $(prefix);\
	fi
	if test ! -d $(bindir); then\
		mkdir $(bindir);\
	fi
	cp $(EXECUTABLES) $(bindir)

.f.o:
	$(F77_COMP) $*.f

.c.o:
	cat $*.c > $*_.c
	$(CC_COMP) -c $*_.c
	mv $*_.o $*.o
	/bin/rm $*_.c


OBJECTS = wave.o init_f.o init_f_t.o u_f.o u_f_t.o res_f.o res_f_t.o
wave: $(OBJECTS)
	$(CC_LOAD) $(OBJECTS) $(LIBS) -o wave

########################################################################
# Clean-up
########################################################################
clean:
	/bin/rm $(EXECUTABLES)
	/bin/rm *_.c > /dev/null 2>&1
	/bin/rm *.o
	/bin/rm *.a
	/bin/rm DISPERSAL BH
	/bin/rm tau maxtmovr

vclean: clean
	/bin/rm  Makefile config.cache config.log config.status
	/bin/rm *.sdf
	/bin/rm -rf autom4te.cache
