CC := g++

STD := -std=c++11

LIB := -lblas -llapack -ltiff

INC := -I ../include/

all:iio.o
	$(CC) iio.o testNystrom.cpp $(INC) $(STD) $(LIB) -o a.out

iio.o:../src/iio.c
	gcc ../src/iio.c -o iio.o -c $(INC) -std=c99

