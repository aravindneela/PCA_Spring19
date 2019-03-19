
CC=mpicc

CFLAGS= -g -Wall -O2


TARGET=cmtbonebe

all: $(TARGET)

$(TARGET): main.o dstructs.o flux.o params.o
	$(CC) -o $@ $^

main.o: main.c dstructs.h utils.h params.h flux.h
	$(CC) -c $(CFLAGS) main.c

flux.o: flux.c flux.h dstructs.h params.h
	$(CC) -c $(CFLAGS) flux.c

dstructs.o: dstructs.c dstructs.h params.h
	$(CC) -c $(CFLAGS) dstructs.c

params.o: params.c params.h
	$(CC) -c $(CFLAGS) params.c

clean:
	rm -rf *.o $(TARGET)
