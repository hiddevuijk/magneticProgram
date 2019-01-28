TARGET = test.exe
OBJS = main.o
CC = g++
OPTFLAG =  -O3 -DNDEBUG
CFLAGS = -c -Wall -Og $(OPTFLAG)
LFLAGS = -Wall -Og

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS)  $(OBJS) -o $(TARGET)

main.o: main.cpp 
	$(CC) $(CFLAGS) main.cpp


.PHONY: clean
clean:
	rm -f  $(OBJS) $(TARGET)

