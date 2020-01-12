T0 = dbAssignment

CC = g++
CFLAGS = -g -std=c++17 -c -Wall -O2
LDFLAGS =

OBJ_T0 = dbAssignment.o

all: $(T0)
# $(T6) $(T3) $(T4)

$(T0): $(OBJ_T0)
	$(CC) $(LDFLAGS) -o $@ $^

.cc.o:
	$(CC) $(CFLAGS) $<

clean:
	rm -f *~ *.o *.exe *.stackdump $(T0)
