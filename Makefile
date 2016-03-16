CFLAGS   = -I./topspek_functions -I./gnuplot_i -I./lin_eq_solver -I./utils -o2 -Wall

all: gnuplot_i.o lin_eq_solver.o dynamic_arrays.o peak_find.o topspek_functions.o topspek

topspek: topspek.c topspek.h peak_find.o topspek_functions.o
	gcc topspek.c $(CFLAGS) -o topspek gnuplot_i.o lin_eq_solver.o dynamic_arrays.o peak_find.o topspek_functions.o
	
gnuplot_i.o: gnuplot_i/gnuplot_i.c gnuplot_i/gnuplot_i.h
	gcc -I./gnuplot_i -o2 -c -o gnuplot_i.o gnuplot_i/gnuplot_i.c
	
lin_eq_solver.o:lin_eq_solver/lin_eq_solver.c lin_eq_solver/lin_eq_solver.h
	gcc -I./lin_eq_solver -o2 -c -o lin_eq_solver.o lin_eq_solver/lin_eq_solver.c
	
dynamic_arrays.o:utils/dynamic_arrays.c utils/dynamic_arrays.h
	gcc -I./utils -o2 -c -o dynamic_arrays.o utils/dynamic_arrays.c

peak_find.o:utils/peak_find.c utils/peak_find.h
	gcc -I./utils -o2 -c -o peak_find.o utils/peak_find.c

topspek_functions.o:topspek_functions/topspek_functions.c topspek_functions/topspek_functions.h
	gcc $(CFLAGS) -c -o topspek_functions.o topspek_functions/topspek_functions.c

clean:
	rm -rf *~ *.o topspek *tmpdatafile*
