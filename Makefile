TARGET = test.exe
OBJS = main.o
CC = g++
OPTFLAG =  -O3 -DNDEBUG
CFLAGS = -c -Wall -Og -g
LFLAGS = -Wall -Og -g
#CFLAGS = -c -Wall -O3 -DNDEBUG
#LFALGS = -Wall -O3 -DNDEBUG

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS)  $(OBJS) -o $(TARGET)

main.o: main.cpp 
	$(CC) $(CFLAGS) main.cpp


.PHONY: clean
clean:
	rm -f  $(OBJS) $(TARGET)

