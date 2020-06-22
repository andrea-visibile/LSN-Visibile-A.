CC = g++
CFLAGS = -O3 --std=c++11

all: esercitazione91.x settings.C esercitazione91.h
	./esercitazione91.x
clean:
	rm -f l2circle.dat
	rm -f l2square.dat
esercitazione91.x: esercitazione91.o random.o settings.C esercitazione91.h
	$(CC) random.o esercitazione91.o -o esercitazione91.x 

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
esercitazione91.o : esercitazione91.C settings.C esercitazione91.h
	$(CC) -c esercitazione91.C -o esercitazione91.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

