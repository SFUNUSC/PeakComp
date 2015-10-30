CFLAGS   = -I./gnuplot_i -o2 -Wall

all: peak_comp gnuplot_i.o

peak_comp: peak_comp.c peak_comp.h read_config.c gnuplot_i.o
	gcc peak_comp.c $(CFLAGS) -o peak_comp gnuplot_i.o
	
gnuplot_i.o: gnuplot_i/gnuplot_i.c gnuplot_i/gnuplot_i.h
	gcc -I./gnuplot_i -o2 -c -o gnuplot_i.o gnuplot_i/gnuplot_i.c

clean:
	rm -rf *~ gnuplot_i.o peak_comp *tmpdatafile*
