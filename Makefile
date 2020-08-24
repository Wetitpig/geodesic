CC=gcc
CFLAGS=-O3 -Iinclude

LDFLAGS=-lm

HEADERS = constants.h functions.h haversine.h vincenty.h

%.o: src/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@

all: src/haversine.o src/vincenty.o src/main.o src/io.o src/math.o
	@mkdir -p bin
	$(CC) $(CFLAGS) -o bin/geodesic $^ $(LDFLAGS)
	@strip bin/geodesic

clean:
	@rm -rf src/*.o
