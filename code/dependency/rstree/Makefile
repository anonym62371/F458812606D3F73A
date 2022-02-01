all: clean all_obj

all_obj: qsort.o poly34.o donut.o rstree_build.o rstree_search.o

qsort.o:
	gcc -c qsort.c -o qsort.o -w

poly34.o:
	gcc -c poly34.c -o poly34.o -w

donut.o:
	gcc -c donut.c -o donut.o -w

rstree_search.o:
	gcc -c rstree_search.c -o rstree_search.o -w

rstree_build.o:
	gcc -c rstree_build.c -o rstree_build.o -w

clean:
	rm -f *.o
	rm -f *.out