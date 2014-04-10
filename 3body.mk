graph.py : data.dat
	python graph.py

data.dat : a.out
	./a.out > data.dat

a.out : main_3body.c
	cc main_3body.c -lm

#clean: 
#	rm data.dat a.out
