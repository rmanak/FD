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
OBJECTS = main.o u_f.o  init_f.o  ire_f.o  res_f.o  res_f_t.o  u_f_t.o  init_f_t.o 
main: $(OBJECTS)
	$(CC_LOAD) $(OBJECTS) $(LIBS) -o main
clean:
	/bin/rm -rf $(EXECUTABLE) $(OBJECTS)
fix:
	sed s/main/main/g < Makefile > Makefile.tmp
	/bin/mv -f Makefile.tmp Makefile
