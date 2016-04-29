.IGNORE:
SHELL = /bin/sh
LIBS = -lm -lbbhutil -lvutil -lifcore 
LDFLAGS = -O3 -L/usr/loca/lib -L.
CC = icc
CFLAGS = -O3
CPPFLAGS = -I. -I/usr/local/include
CC_COMP = $(CC) -c $(CFLAGS) $(CPPFLAGS)
CC_LOAD = $(CC) $(LDFLAGS)
F77 = ifort
F77FLAGS = -O3 -w -w90 -w95 -cm -Vaxlib
F77_COMP = $(F77) -c $(F77FLAGS)
F77_LOAD = $(F77) $(F77FLAGS) $(F77_LDFLAGS) $(LDFLAGS)
EXECUTABLE = main
all: $(EXECUTABLE)
.f.o:
	$(F77_COMP) $*.f
.c.o:
	$(CC_COMP) -c $*.c
OBJECTS = main.o init_pi.o  init_pp.o  ire_pi.o  ire_pp.o  res_pi.o  res_pp.o  u_pi.o  u_pp.o  calc_trr.o  calc_ttt.o  init_phi.o ode_driver.o ire_a.o ire_alpha.o
main: $(OBJECTS)
	$(CC_LOAD) $(OBJECTS) $(LIBS) -o main
clean:
	/bin/rm -rf $(EXECUTABLE) $(OBJECTS)
fix:
	sed s/main/main/g < Makefile > Makefile.tmp
	/bin/mv -f Makefile.tmp Makefile
