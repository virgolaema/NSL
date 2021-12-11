CC = mpic++
CFLAGS = -Wall # -O3  --std=c++11

fast: compila esegui

compila: main.o lib.o random.o TSP.o
	$(CC) main.o lib.o random.o TSP.o -o main
main.o : main.C ../lib.h ../random.h TSP.h
	$(CC) -c main.C 
lib.o : ../lib.C ../lib.h 
	$(CC) -c ../lib.C 
TSP.o : TSP.C TSP.h 
	$(CC) -c TSP.C 
random.o: ../random.cpp ../random.h 
	$(CC) -c ../random.cpp 

esegui: main
	mpirun -np 4 ./main

clean: 
	rm *.o main