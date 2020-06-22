CC = g++
CFLAGS =  -O3 --std=c++11

all: MolDyn_NVE.x 
	./MolDyn_NVE.x

MolDyn_NVE.x: MolDyn_NVE.o random.o
	$(CC) random.o MolDyn_NVE.o -o MolDyn_NVE.x 

MolDyn_NVE.o : MolDyn_NVE.cpp input.dat
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o old.0 old.final config.final
Solid:
	cp inputs/InputArgonSolid.dat input.dat
Liquid:
	cp inputs/InputArgonLiquid.dat input.dat
Gas:
	cp inputs/InputArgonGas.dat input.dat	
generalEq: 
	./clean.sh
	cp inputs/eq0.dat input.dat
	./MolDyn_NVE.x
	cp inputs/eq1.dat input.dat
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	
	
ArgonSolidEq: 
	./clean.sh
	rm -rf Argon/Solid/* 
	cp inputs/inputArgonSolid0.dat input.dat
	./MolDyn_NVE.x
	cp inputs/inputArgonSolid1.dat input.dat
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./clean.sh
	rm -rf Argon/Solid/*
	cp inputs/inputArgonSolid.dat input.dat
	./MolDyn_NVE.x
ArgonLiquidEq:
	rm -rf Argon/Liquid/* 
	./clean.sh
	cp inputs/inputArgonLiquid0.dat input.dat
	./MolDyn_NVE.x
	cp inputs/inputArgonLiquid1.dat input.dat
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./clean.sh
	rm -rf Argon/Liquid/*
	cp inputs/inputArgonLiquid.dat input.dat
	./MolDyn_NVE.x
ArgonGasEq: 
	rm -rf Argon/Gas/*
	./clean.sh
	cp inputs/inputArgonGas0.dat input.dat
	./MolDyn_NVE.x
	cp inputs/inputArgonGas1.dat input.dat
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	./MolDyn_NVE.x
	rm -rf Argon/Gas/*
	./clean.sh
	cp inputs/inputArgonGas.dat input.dat
	./MolDyn_NVE.x
