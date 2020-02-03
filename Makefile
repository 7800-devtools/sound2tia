# The sound2tia generic makefile

CC     = gcc
LIBS   = -lsndfile -lfftw3
DEBUG  = -Wall -g
 
all: sound2tia
 
sound2tia:	sound2tia.c tiavals.h
	$(CC) sound2tia.c $(DEBUG) $(LIBS) -o sound2tia

clean:
	rm -f *.o sound2tia
