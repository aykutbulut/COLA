CC=g++
PROGNAME=cut
INCLUDES=`pkg-config --cflags cola` -I/home/aykut/research/OsiConic/build/include
CFLAGS=-g
LDFLAGS=`pkg-config --libs cola`
LIBS =

OBJECTS=cut.o

all: $(PROGNAME)

$(PROGNAME): $(OBJECTS)
	$(CC) -o $(PROGNAME) $(OBJECTS)  $(LDFLAGS) $(LIBS)

$(OBJECTS): Makefile

.cpp.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $@ $<

clean:
	rm *.o $(PROGNAME)

test: $(PROGNAME)
	./$(PROGNAME) mir_example.mps
