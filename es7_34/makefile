CC = g++
CFLAGS = -Wall -O3 --std=c++11

fast: compila esegui

compila : MolDyn_NVE.o 
	$(CC)  MolDyn_NVE.o -o MolDyn_NVE.exe
main.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm *.o *.out MolDyn_NVE.exe
esegui:
	./MolDyn_NVE.exe