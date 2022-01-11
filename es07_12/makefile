CC = g++
CFLAGS = -Wall -O3 --std=c++11

fast: pulisci compila esegui

compila: Monte_Carlo_NVT.o random.o
	$(CC) random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

esegui: Monte_Carlo_NVT.exe
	./Monte_Carlo_NVT.exe

clean :
	rm *.o Monte_Carlo_NVT.exe

liquid:
	cp input.liquid input.dat
	cp config.liquid config.0

solid:
	cp input.solid input.dat
	cp config.solid config.0

gas: 
	cp input.gas input.dat
	cp config.gas config.0

fcc: 
	cp config.fcc config.dat

pulisci: 
	rm output*