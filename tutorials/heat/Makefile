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
OBJECTS = init_f.o    update_f.o main.o
FORTS = init_f.f  update_f.f
CALLS = init_f_call   update_f_call
HEADS = init_f.h  update_f.h

all: $(FORTS) $(OBJECTS) $(EXECUTABLE)

init_f.f: heat.mpl
	maple < heat.mpl > /dev/null 2>&1
	sleep 1
	ls -ltr *.f *.h *_call
.f.o:
	$(F77_COMP) $*.f
.c.o:
	$(CC_COMP)  $*.c

main: $(OBJECTS)
	$(CC_LOAD) $(OBJECTS) $(LIBS) -o main

clean:
	/bin/rm $(EXECUTABLE) $(OBJECTS) 
vclean:
	/bin/rm -rf *~ *_ *.sdf $(EXECUTABLE) $(OBJECTS) $(FORTS) $(CALLS) $(HEADS)
fix:
	sed s/main/auto/g < Makefile > Makefile.tmp
	/bin/mv -f Makefile.tmp Makefile
