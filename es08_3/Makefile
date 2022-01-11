LIBS:=`root-config --libs`
INCS:=`root-config --cflags`


%.o : %.cpp
	g++ -Wall -c $< ${INCS}

qmc1d: qmc1d.o
	g++ -O3 -Wall -o $@ $^ ${LIBS}

esegui: qmc1d
	./qmc1d

clean:
	rm *.o qmc1d potential.dat kinetic.dat probability.dat

pigs:
	cp input.pigs input.dat

pimc:
	cp input.pimc input.dat

