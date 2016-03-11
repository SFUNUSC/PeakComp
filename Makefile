CFLAGS   = -I./peakcomp_functions -I./gnuplot_i -I./lin_eq_solver -I./utils -o2 -Wall

all: gnuplot_i.o lin_eq_solver.o dynamic_arrays.o peak_find.o peakcomp_functions.o topspek

topspek: peak_comp.c peak_comp.h peak_find.o peakcomp_functions.o
	gcc peak_comp.c $(CFLAGS) -o topspek gnuplot_i.o lin_eq_solver.o dynamic_arrays.o peak_find.o peakcomp_functions.o
	
gnuplot_i.o: gnuplot_i/gnuplot_i.c gnuplot_i/gnuplot_i.h
	gcc -I./gnuplot_i -o2 -c -o gnuplot_i.o gnuplot_i/gnuplot_i.c
	
lin_eq_solver.o:lin_eq_solver/lin_eq_solver.c lin_eq_solver/lin_eq_solver.h
	gcc -I./lin_eq_solver -o2 -c -o lin_eq_solver.o lin_eq_solver/lin_eq_solver.c
	
dynamic_arrays.o:utils/dynamic_arrays.c utils/dynamic_arrays.h
	gcc -I./utils -o2 -c -o dynamic_arrays.o utils/dynamic_arrays.c

peak_find.o:utils/peak_find.c utils/peak_find.h
	gcc -I./utils -o2 -c -o peak_find.o utils/peak_find.c

peakcomp_functions.o:peakcomp_functions/peakcomp_functions.c peakcomp_functions/peakcomp_functions.h
	gcc $(CFLAGS) -c -o peakcomp_functions.o peakcomp_functions/peakcomp_functions.c

clean:
	rm -rf *~ *.o topspek *tmpdatafile*
