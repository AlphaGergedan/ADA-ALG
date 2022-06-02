.PHONY: all
.PHONY: clean

all: tests

quicksort.o:
	g++ ../src/quicksort.cpp -c -o ../build/quicksort.o

tests: quicksort.o
	g++ -Wall -g -pthread quicksort-tests.cpp ../src/quicksort.cpp -lgtest_main  -lgtest -lpthread -o tests

clean:
	rm tests
