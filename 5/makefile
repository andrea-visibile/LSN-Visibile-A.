CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercitazione51.x 
	./esercitazione51.x
	


esercitazione51.x: esercitazione51.o random.o 
	$(CC) random.o esercitazione51.o -o esercitazione51.x 

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
esercitazione51.o : esercitazione51.C esercitazione51.h settings.C
	$(CC) -c esercitazione51.C -o esercitazione51.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
