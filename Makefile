CC = gcc

main: main.c
	$(CC) -o main main.c -fopenmp

parallel: parallel.c
	$(CC) -o parallel parallel.c -fopenmp

clean:
	-rm -f main.o
	-rm -f main
	-rm -f parallel.o
	-rm -f parallel