TARGET = test.exe
OBJS = main.o
CC = g++
OPTFLAG =  -O3 -DNDEBUG
CFLAGS = -c -Wall  -g
LFLAGS = -Wall  -g
#CFLAGS = -c -Wall -O3 -DNDEBUG
#LFALGS = -Wall -O3 -DNDEBUG

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS)  $(OBJS) -o $(TARGET)

main.o: main.cpp xyz.h bfield.h system_AT.h\
	 system.h density.h orientation.h flux.h
	$(CC) $(CFLAGS) main.cpp


.PHONY: clean
clean:
	rm -f  $(OBJS) $(TARGET)

