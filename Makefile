OBJS = parseinput.o kdtree.o kdtree_greedy.o kdtree_localsearch.o \
       kdtree_semigreedy.o util.o dist.o heapmax.o naive.o naive_grasp.o
TARGET = tsp test_kdtree
LIBS = -lm
MYFLAGS = -Wall
CFLAGS =
CC = gcc

all: $(TARGET)

tsp: $(OBJS) main.o
	$(CC) $(MYFLAGS) $(CFLAGS) $(LIBS) -o $@ $(OBJS) main.c

test_kdtree: $(OBJS) test_kdtree.o
	$(CC) $(MYFLAGS) $(CFLAGS) $(LIBS) -o $@ $(OBJS) test_kdtree.c

%.o: %.c
	$(CC) $(MYFLAGS) $(CFLAGS) -c $<

clean:
	rm -rf $(OBJS) $(TARGET) main.o test_kdtree.o
