.PHONY: all
.PHONY: clean

all: extended_euclid

extended_euclid.o:
	g++ src/extended_euclid.cpp -c -o build/extended_euclid.o

extended_euclid: extended_euclid.o
	g++ src/main.cpp build/extended_euclid.o -o build/extended_euclid

clean:
	rm build/*
