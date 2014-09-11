CC=g++ -fopenmp -Wall -O4
OBJS=meteor.o
LFLAGS= -lm -L /home/ram80/include -lgomp
CFLAGS= -I /home/ram80/include -lgomp
all: meteor
meteor: $(OBJS)
	$(CC) -o meteor $(OBJS) $(CFLAGS)
meteor.o: meteor-fork.cpp
	$(CC) -o meteor.o -c meteor-fork.cpp $(LFLAGS)
clean:
	rm *.o meteor
