CC = g++
CFLAGS = -O3 --std=c++11

all: esercitazione101.x 
	./esercitazione101.x


esercitazione101.x: esercitazione101.o random.o
	$(CC) random.o esercitazione101.o -o esercitazione101.x 
 

esercitazione101.o : esercitazione101.C esercitazione101.h
	$(CC) -c esercitazione101.C -o esercitazione101.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
