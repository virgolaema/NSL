CC = g++
CFLAGS = -Wall # -O3  --std=c++11

MolDyn_NVE:  lib.o random.o MolDyn_NVE.o  
	$(CC) MolDyn_NVE.o lib.o random.o -o MolDyn_NVE
prepare: prepare.sh
	./prepare.sh
config.0: #copying the config to start with
	cp config.final config.0
MolDyn_NVE.o : MolDyn_NVE.C ../lib.h ../random.h MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.C 
lib.o : ../lib.C ../lib.h 
	$(CC) -c ../lib.C 
random.o: ../random.cpp ../random.h 
	$(CC) -c ../random.cpp 

esegui: MolDyn_NVE 
	./MolDyn_NVE

clean: 
	rm *.o MolDyn_NVE

liquid: 
	cp inputs/input.liquid input.dat

solid: 
	cp inputs/input.solid input.dat

gas: 
	cp inputs/input.gas input.dat

fcc: 
	cp configs/config.fcc config.0