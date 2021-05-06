CC = g++
CFLAGS = -Wall # -O3  --std=c++11

fast: compila esegui

compila: main.o lib.o random.o 
	$(CC) main.o lib.o random.o -o main
main.o : main.C ../lib.h ../random.h
	$(CC) -c main.C 
lib.o : ../lib.C ../lib.h 
	$(CC) -c ../lib.C 
random.o: ../random.cpp ../random.h 
	$(CC) -c ../random.cpp 

esegui: main
	./main

clean: 
	rm *.o main

U1: 
	cp U1.dat input.dat
G1: 
	cp G1.dat input.dat
U2: 
	cp U2.dat input.dat
G2: 
	cp G2.dat input.dat
