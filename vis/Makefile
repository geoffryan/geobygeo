
CC = gcc

CFLAGS ?= -Wall
CURL ?= curl
UNTAR ?= tar -xvf
CD ?= cd
RM ?= rm -f
OS ?= MacOSX10.6.8

LIBS = -framework OpenGL -framework GLUT
INCS = -I/usr/X11R6/include

default: bin/geovis

bin/geovis: geovis.c
	@mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $^ $(INCS) $(LIBS)
	
clean:
	$(RM) bin/geovis
	
.FORCE:

