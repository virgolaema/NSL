CC = g++
CFLAGS = -Wall -O3 --std=c++11

fast: pulisci prec compila esegui

test: prec compila esegui 

compila: main.o random.o lib.o
	$(CC) random.o lib.o main.o -o main
main.o : main.C ../lib.h ../random.h
	$(CC) -c main.C 
lib.o : ../lib.C ../lib.h 
	$(CC) -c ../lib.C 
random.o: ../random.cpp ../random.h 
	$(CC) -c ../random.cpp 

esegui: main
	./main

prec: 
	cp config.final config.0

pulisci: 
	rm *out*

clean:
	rm *.o main

sposta: 
	mv *.out results
M:
	cp inputs/metropolis.dat input.dat

G:
	cp configs/gibbs.dat input.dat

M_m:
	cp configs/metropolis_m.dat input.dat

G_m:
	cp configs/gibbs_m.dat input.dat