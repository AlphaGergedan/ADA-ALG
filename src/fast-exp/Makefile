.PHONY: all
.PHONY: clean

all: fast_exp

fast_exp.o:
	g++ src/fast_exp.cpp -c -o build/fast_exp.o

fast_exp: fast_exp.o
	g++ src/main.cpp build/fast_exp.o -o build/fast_exp

clean:
	rm build/*




