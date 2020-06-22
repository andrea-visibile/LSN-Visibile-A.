CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o random.o
	$(CC) random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_NVT.exe

input_gas:
	cp input.gas input.dat
	./Monte_Carlo_NVT.exe
	mv output.* gas
	mv pot_ist.dat gas
	mv press_ist.dat gas

input_liquid:
	cp input.liquid input.dat
	./Monte_Carlo_NVT.exe
	mv output.* liquid
	mv pot_ist.dat liquid
	mv press_ist.dat liquid
input_solid:
	cp input.solid input.dat
	./Monte_Carlo_NVT.exe
	mv output.* solid	
	mv pot_ist.dat solid
	mv press_ist.dat solid
