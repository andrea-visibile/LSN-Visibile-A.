CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercitazione21.x esercitazione22.x 
	./esercitazione21.x
	./esercitazione22.x


esercitazione21.x: esercitazione21.o random.o
	$(CC) random.o esercitazione21.o -o esercitazione21.x 
esercitazione22.x: esercitazione22.o random.o
	$(CC) random.o esercitazione22.o -o esercitazione22.x 
main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
esercitazione21.o : esercitazione21.C esercitazione21.h
	$(CC) -c esercitazione21.C -o esercitazione21.o $(CFLAGS)
esercitazione22.o : esercitazione22.C esercitazione22.h
	$(CC) -c esercitazione22.C -o esercitazione22.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
