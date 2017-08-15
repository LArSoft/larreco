# Macros
CC = g++
CCflags = -Wall

all: tc_sh_debugger

clean:
	rm -f tc_sh_debugger *.o *~

tc_sh_debugger: showertree.C showertree.h
	$(CC) $(CCFlags) -o tc_sh_debugger showertree.C `root-config --cflags --glibs`