graph.py : data.dat
	python graph.py

data.dat : a.out main_3body.c
	a.out > data.dat
	cc main_3body.c -lm


