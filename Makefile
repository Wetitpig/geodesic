CC=gcc
CFLAGS=-O3 -Iinclude

LDFLAGS=-lm

all: haversine vincenty

haversine:
	@mkdir -p bin
	$(CC) $(CFLAGS) -o bin/haversine src/haversine.c src/math.c src/io.c $(LDFLAGS)

vincenty:
	@mkdir -p bin
	$(CC) $(CFLAGS) -o bin/vincenty src/vincenty.c src/math.c src/io.c $(LDFLAGS)

clean:
	rm -rf bin
