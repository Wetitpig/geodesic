CC=gcc
CFLAGS=-O3

LDFLAGS=-lm

all: haversine vincenty

haversine:
	$(CC) $(CFLAGS) -o haversine haversine.c math.c io.c $(LDFLAGS)

vincenty:
	$(CC) $(CFLAGS) -o vincenty vincenty.c math.c io.c $(LDFLAGS)

clean:
	rm -f haversine vincenty
