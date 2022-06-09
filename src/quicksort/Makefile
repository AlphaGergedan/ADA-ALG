.PHONY: all
.PHONY: clean

all: quicksort

quicksort.o:
	g++ src/quicksort.cpp -c -o build/quicksort.o

quicksort: quicksort.o
	g++ src/main.cpp build/quicksort.o -o build/quicksort

clean:
	rm build/*
