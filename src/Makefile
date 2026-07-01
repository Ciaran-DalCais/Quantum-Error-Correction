CC = gcc

QUEST_ROOT = /home/gandalf/Code/Github/QuEST
QUEST_BUILD = $(QUEST_ROOT)/build

CFLAGS = -Wall -g \
          -I$(QUEST_ROOT)/quest/include \
          -I$(QUEST_ROOT) \
          -I$(QUEST_BUILD)

LDFLAGS = -L$(QUEST_BUILD) \
          -Wl,-rpath,$(QUEST_BUILD)

LDLIBS = -lQuEST -lm -llapacke -llapack -lblas 

SRC = main.c encoder_functions.c 
OBJ = $(SRC:.c=.o)

main: $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) $(LDLIBS) -o main

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o main

.PHONY: all clean
