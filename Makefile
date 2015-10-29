CFLAGS   = -O -Wall

all: peak_comp

peak_comp: peak_comp.c
	gcc peak_comp.c $(CFLAGS) -o peak_comp

clean:
	rm -rf *~ peak_comp
