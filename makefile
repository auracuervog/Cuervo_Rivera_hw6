all:	Presas.pdf

Presas.pdf:	graph.py
	python graph.py
graph.py:	presa.dat

presa.dat:	bla.out
	./bla.out 29 20
bla.out:	punto_1.c
	cc punto_1.c -o bla.out
