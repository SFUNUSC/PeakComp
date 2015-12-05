CFLAGS   = -I./gnuplot_i -I./lin_eq_solver -o2 -Wall

all: peak_comp gnuplot_i.o lin_eq_solver.o

peak_comp: peak_comp.c peak_comp.h read_config.c gnuplot_i.o lin_eq_solver.o
	gcc peak_comp.c $(CFLAGS) -o peak_comp gnuplot_i.o lin_eq_solver.o
	
gnuplot_i.o: gnuplot_i/gnuplot_i.c gnuplot_i/gnuplot_i.h
	gcc -I./gnuplot_i -o2 -c -o gnuplot_i.o gnuplot_i/gnuplot_i.c
	
lin_eq_solver.o:lin_eq_solver/lin_eq_solver.c lin_eq_solver/lin_eq_solver.h
	gcc -I./lin_eq_solver -o2 -c -o lin_eq_solver.o lin_eq_solver/lin_eq_solver.c

clean:
	rm -rf *~ *.o peak_comp *tmpdatafile*
