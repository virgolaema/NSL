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
