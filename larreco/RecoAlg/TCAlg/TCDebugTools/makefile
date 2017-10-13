# Macros
CC = g++
CCflags = -Wall

all: rsf_evd

clean:
	rm -f rsf_evd *.o *~

rsf_evd: showertree.C showertree.h
	$(CC) $(CCFlags) -o rsf_evd showertree.C `root-config --cflags --glibs`