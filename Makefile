CC=gcc
CFLAGS=-O3 -Iinclude -g

LDFLAGS=-lm

HEADERS = constants.h functions.h haversine.h vincenty.h greatcircle.h karney.h karney_lookup.h

%.o: src/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@

all: debug
	@strip bin/geodesic

debug: src/haversine.o src/vincenty.o src/main.o src/io.o src/math.o src/greatcircle.o src/karney.o
	@mkdir -p bin
	$(CC) $(CFLAGS) -o bin/geodesic $^ $(LDFLAGS)

install: all
	@mv bin/geodesic $(prefix)/bin/geodesic
	rm -rf bin

clean:
	@rm -rf src/*.o
