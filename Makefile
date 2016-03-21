CFLAGS   = -I./topspek_functions -I./gnuplot_i -I./lin_eq_solver -I./utils -o2 -Wall

all: lib topspek

topspek: topspek.c topspek.h peak_find.o
	@echo Making topspek...
	gcc topspek.c $(CFLAGS) -o topspek gnuplot_i.o lin_eq_solver.o dynamic_arrays.o peak_find.o
	@echo Tidying up...
	rm -rf *~ *.o
	
lib: gnuplot_i/gnuplot_i.c gnuplot_i/gnuplot_i.h lin_eq_solver/lin_eq_solver.c lin_eq_solver/lin_eq_solver.h utils/dynamic_arrays.c utils/dynamic_arrays.h utils/peak_find.c utils/peak_find.h
	@echo Making libraries...
	gcc -I./gnuplot_i -o2 -c -o gnuplot_i.o gnuplot_i/gnuplot_i.c
	gcc -I./lin_eq_solver -o2 -c -o lin_eq_solver.o lin_eq_solver/lin_eq_solver.c
	gcc -I./utils -o2 -c -o dynamic_arrays.o utils/dynamic_arrays.c
	gcc -I./utils -o2 -c -o peak_find.o utils/peak_find.c

clean:
	@echo Cleaning up...
	rm -rf *~ *.o topspek *tmpdatafile*
