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
OBJECTS = main.o u_f.o  init_f.o  u_f_t.o  init_f_t.o res_f.o res_f_t.o ire_f.o
FORTS= u_f.f init_f.f u_f_t.f init_f_t.f res_f.f res_f_t.f ire_f.f

all: init_f.f $(EXECUTABLE)

init_f.f: w1dpb.mpl
	maple < w1dpb.mpl > /dev/null 2>&1
	sleep 2
	ls -ltr *.f *.h *_call

.f.o:
	$(F77_COMP) $*.f

.c.o:
	$(CC_COMP) -c $*.c

main: $(OBJECTS) init_f.f
	$(CC_LOAD) $(OBJECTS) $(LIBS) -o main
clean:
	/bin/rm -rf $(EXECUTABLE) $(OBJECTS)
vclean: clean
	/bin/rm -rf *~ *_ *.sdf $(FORTS) 
fix:
	sed s/main/main/g < Makefile > Makefile.tmp
	/bin/mv -f Makefile.tmp Makefile
