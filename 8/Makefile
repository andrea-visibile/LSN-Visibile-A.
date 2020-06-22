CC = g++
CFLAGS = -O3 --std=c++11

all: esercitazione81.x
	rm -rf H.dat
	rm -rf x.dat
	./esercitazione81.x


esercitazione81.x: esercitazione81.o random.o
	$(CC) random.o esercitazione81.o -o esercitazione81.x 

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
esercitazione81.o : esercitazione81.C 
	$(CC) -c esercitazione81.C -o esercitazione81.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
