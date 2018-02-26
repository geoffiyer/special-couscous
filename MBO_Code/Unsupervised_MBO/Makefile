CC := g++

SRCDIR := src

BINDIR := bin

BUILDDIR := build

TARGET := $(BINDIR)/b.exe


SRCEXT := cpp

SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))

OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o)) iio.o

CFLAGS := -O3

LIB := -lblas -llapack -ltiff

OMPOPTS := -fopenmp -DOMP

INC := -I include

$(TARGET):$(OBJECTS)
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(OMPOPTS) $^ $(LIB) -o $(TARGET) 

$(BUILDDIR)/%.o:$(SRCDIR)/%.$(SRCEXT) 
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) $(OMPOPTS) $(LIB) $(INC) -c -o $@ $^
iio.o:$(SRCDIR)/iio.c
	gcc $(SRCDIR)/iio.c -o iio.o -c $(INC) -std=c99
clean:
	$(RM) -r $(BUILDDIR) $(BINDIR) iio.o classification_result.txt
.PHONY: clean
