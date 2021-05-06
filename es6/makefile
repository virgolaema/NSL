CC = g++
CFLAGS = -Wall -O3 --std=c++11

fast: pulisci compila esegui

compila: main.o random.o lib.o
	$(CC) random.o lib.o main.o -o main
main.o : main.C main.h ../lib.h ../random.h
	$(CC) -c main.C 
lib.o : ../lib.C ../lib.h 
	$(CC) -c ../lib.C 
random.o: ../random.cpp ../random.h 
	$(CC) -c ../random.cpp 

esegui: main
	./main

pulisci: 
	rm *.0

clean:
	rm *.o main
