#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*
*Metodo principal del programa, llama a los otros metodos y genera cinco archivos de texto con los datos de la evolucion de las particulas
*Entradas: -argc: tamanio del vector de parametros
*-*argv[]: es el vector de parametros, con el archivo de posicions generado en el primer programa y el tiempo T en que se deja evolucionar el sistema
*Salidas: -0
*/
double *Evolucion (double pos0x, double pos0y, double vel0x, double vel0y);

int main (int argc, char *argv[])
{
	FILE *Aentrada;
	FILE *a1;
	FILE *a2;
	FILE *a3;
	FILE *a4;
	FILE *a;
	Resolver(0.0,0.0,0.0,0.0);

	return 0;
}

/*
*Metodo para resolver la ecuacion diferencial de segundo orden para una particula orbitando el centro de la galaxia. Utiliza el codigo de Runge-Kutta de cuarto orden disponible en el repositorio
*Entradas: Condiciones iniciales: pos0x: posicion inicial en x; pos0y: posicion inicial en y; vel0x: velocidad inicial en x; viy: velocidad inicial en y; tiempo T en el que se va a evaluar la evolucion del sistema.
*Salidas: *SolPos: vector con las posiciones de las particulas en funcion del tiempo
*-*SolVel: vector con las velocidades de las particulas en funcion del tiempo
*/

double *Movimiento (double pos0x, double pos0y, double vel0x, double vel0y, int T)
{
	double *SolPos;
	double *SolVel;
	SolPos = malloc(15 * sizeof(int));
	SolVel = malloc(15 * sizeof(int));
	return Pos_ev, Vel_ev;

}
