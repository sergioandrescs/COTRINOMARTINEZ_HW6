all:  depredador_presa.o datos.o graficas.o

graficas.o: graficas.py
	python graficas.py 

datos.o: volterra_lotka.x
	./volterra_lotka.x

depredador_presa.o: volterra_lotka.c
	cc volterra_lotka.c -o volterra_lotka.x

clean:
	rm *.x *.txt *.png
