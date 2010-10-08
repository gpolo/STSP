OBJS = parseinput.o kdtree.o kdtree_greedy.o kdtree_localsearch.o \
       kdtree_semigreedy.o util.o dist.o heapmax.o naive.o naive_grasp.o
TARGET = tsp test_kdtree
LIBS = -lm
MYFLAGS =
CFLAGS =
CC = gcc

all: $(TARGET)

tsp: $(OBJS) main.o
	$(CC) -o $@ $(OBJS) main.c $(MYFLAGS) $(CFLAGS) $(LIBS)

test_kdtree: $(OBJS) test_kdtree.o
	$(CC) -o $@ $(OBJS) test_kdtree.c $(MYFLAGS) $(CFLAGS) $(LIBS)

%.o: %.c
	$(CC) -c $< $(MYFLAGS) $(CFLAGS)

clean:
	rm -rf $(OBJS) $(TARGET) main.o test_kdtree.o

clean-profile:
	rm -rf *.gcda *.gcno
