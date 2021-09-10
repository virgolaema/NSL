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

liquid: l_input l_config

l_input:
	cp input.liquid input.dat

l_config:
	cp config.liquid config.0

solid: s_input s_config

s_input:
	cp input.solid input.dat

s_config:
	cp config.solid config.0

gas: g_input g_config

g_input:
	cp input.gas input.dat

g_config:
	cp config.gas config.0

fcc: 
	cp config.fcc config.dat

pulisci: 
	rm output*