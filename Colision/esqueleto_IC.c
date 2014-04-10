#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double *pos(int r);
double *vel(int r);
/*
*Metodo principal del programa que genera un archivo con los datos de posiciones y velocidades de las particulas de la galaxia a partir del uso de pos y vel
*Entradas: argc: tamanio del vector de parametros; *argv[]: es el vector de parametros, con las posiciones y velocidades inicales de la masa central 
*Salidas: archivo de texto con las condiciones iniciales
*/
int main (int argc, char *argv[])
{
	FILE *archivo;
	int R;
	double *Posi;
	double *Velo;
	Posi = pos(R);
	Velo = vel(R);
	return 0;
}
/*
*Metodo para calcular las posiciones iniciales de las particulas
*Entradas: R: entero que indica el radio en kiloparsecs de la galaxia; N: numero de particulas
*Salidas: *pos1: vector con las coordenadas rectangulares de las particulas. La masa central es el origen y las particulas se distribuyen uniformemente de manera aleatoria dentro del circulo de radio R 
*/
double *pos (int R, int N)
{
	double *pos1;
	int N;
	pos1=10
	return pos1;
}

double *vel (double pos1)
{
/*
*Metodo para calcular las velocidades iniciales d elas particulas
*Entradas: pos1
*Salidas: *Velocidades: vector con las componentes rectangulares de la velocidad de cada particula. Se calculan considerando que deben estar en equilibrio, por lo que la aceleracion centripeta se iguala a la fuerza radial ejercidad por la gravedad de la masa central.
*/
	double *vel1;
	vel1=10
	return vel1;

}
