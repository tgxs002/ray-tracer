CC = /usr/bin/g++

#COMPILER FLAGS
CCFLAGS := -Ofast -fopenmp -std=c++11

#include directories
#should include gl.h glut.h etc...
INCDIR := -I/usr/include 
LDLIBS := $(GLLIB)

TARGET = main
OBJS = main.o tracer.o elements.o sampler.o


all: $(TARGET)


$(TARGET): $(OBJS)
	$(CC)  $^ $(CCFLAGS) $(LDLIBS)  -o $@

%.o : %.cpp
	$(CC) $(CCFLAGS) -o $@ -c $(LDLIBS) $(INCDIR) $<

clean:
	rm -f $(OBJS) $(TARGET)

