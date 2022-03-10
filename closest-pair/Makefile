.PHONY: all
.PHONY: closest-pair
.PHONY: main
.PHONY: test
.PHONY: clean

FLAGS=-O0 -g

all: closest-pair main

point.o:
	g++ ${FLAGS} src/point2D.cpp -c -o build/point2D.o

closest-pair.o:
	g++ ${FLAGS} src/closestPair.cpp -c -o build/closestPair.o

closest-pair: point.o closest-pair.o
	g++ ${FLAGS} build/point2D.o build/closestPair.o -o build/closestPair

main:

test: point.o
	g++ ${FLAGS} src/test.cpp build/point2D.o -o build/test

clean:
	rm build/*
