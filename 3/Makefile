CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercitazione31.x 
	./esercitazione31.x
	


esercitazione31.x: esercitazione31.o random.o
	$(CC) random.o esercitazione31.o -o esercitazione31.x 

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
esercitazione31.o : esercitazione31.C esercitazione31.h
	$(CC) -c esercitazione31.C -o esercitazione31.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
