CC=g++ -std=c++11 -O3 -fopenmp

all: FASI LPGO QueryGenerator LabelGenerator

FASI.o: FASI.cpp
	$(CC) -c -I./ FASI.cpp -o FASI.o

LPGO.o: LPGO.cpp
	$(CC) -c -I./ LPGO.cpp -o LPGO.o

LabelGenerator.o: LabelGenerator.cpp
	$(CC) -c -I./ LabelGenerator.cpp -o LabelGenerator.o

QueryGenerator.o: QueryGenerator.cpp
	$(CC) -c -I./ QueryGenerator.cpp -o QueryGenerator.o

LPCSR.o: LPCSR.cpp LPCSR.h Graph.o
	$(CC) -c LPCSR.cpp -o LPCSR.o 

Graph.o: Graph.cpp Graph.h
	$(CC) -c Graph.cpp -o Graph.o

FASI: FASI.o LPCSR.o Graph.o
	$(CC) FASI.o LPCSR.o Graph.o -o FASI

LPGO: LPGO.o LPCSR.o Graph.o
	$(CC) LPGO.o LPCSR.o Graph.o -o LPGO

QueryGenerator: QueryGenerator.o
	$(CC) QueryGenerator.o -o QueryGenerator

LabelGenerator: LabelGenerator.o
	$(CC) LabelGenerator.o -o LabelGenerator

clean:
	rm -rf *.o FASI LPGO QueryGenerator LabelGenerator
