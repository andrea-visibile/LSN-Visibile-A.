CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: esercitazione11.x esercitazione12.x esercitazione13.x
	./esercitazione11.x
	./esercitazione12.x
	./esercitazione13.x

esercitazione11.x: esercitazione11.o random.o
	$(CC) random.o esercitazione11.o -o esercitazione11.x 
esercitazione12.x: esercitazione12.o random.o
	$(CC) random.o esercitazione12.o -o esercitazione12.x 
esercitazione13.x: esercitazione13.o random.o
	$(CC) random.o esercitazione13.o -o esercitazione13.x 
main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
esercitazione11.o : esercitazione11.C esercitazione11.h
	$(CC) -c esercitazione11.C -o esercitazione11.o $(CFLAGS)
esercitazione12.o : esercitazione12.C esercitazione12.h
	$(CC) -c esercitazione12.C -o esercitazione12.o $(CFLAGS)
esercitazione13.o : esercitazione13.C esercitazione13.h
	$(CC) -c esercitazione13.C -o esercitazione13.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
