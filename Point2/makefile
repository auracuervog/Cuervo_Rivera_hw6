all: trayectoria_100000000.0_30.0_2d.pdf trayectoria_100000000.0_30.0_3d.pdf

trayectoria_100000000.0_30.0_2d.pdf:  trayectoria_100000000.0_30.0.dat
	python graph.py  trayectoria_100000000.0_30.0.dat
trayectoria_100000000.0_30.0_3d.pdf:  trayectoria_100000000.0_30.0.dat
	python graesteb.py  trayectoria_100000000.0_30.0.dat
trayectoria_100000000.0_30.0.dat: a.out
	./a.out 100000000 30

a.out: punto2.c
	cc punto2.c -lm
