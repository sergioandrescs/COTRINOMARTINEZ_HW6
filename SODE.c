#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MSOLAR 1.989E30
#define UA 1.496E11
#define YEAR_SEC 31536000
#define G 6.67E-11


/*Solucion de ecuaciones diferenciales con Runge Kutta cuarto orden*/

int main (){
  /*Define todas las variables a usar
    proceso entre x_min y x_max
    k1, k2, k3, k4, k son las pendientes para aplicar este metodo
    x_op y y_op son variable usadas para calcular las pendientes
  */

  double h;
  double y0;
  double x_min;
  double x_max;
  double k1, k2, k3, k4, k;
  int n_data,i;
  double *d;  
  double x_op, y_op; 
  double new_G =(G*MSOLAR*pow(YEAR_SEC,2))/pow(UA,3)

  x_max=1.5;
  x_min=0.0;
  h=0.01;
  n_data=((x_max-x_min)/h);

  d=malloc(n_data*2*sizeof(double));

  d[0]=x_min;
  d[n_data]=0;

  double fun1 (double x, double y1, double y2){
    return ;
  }
  double fun2 (double x, double x1, double y2 ){
    return 5;
  }


  for (i=1;i<n_data;i++){
    /*Primera pendiente*/   
    k1=func(d[0],d[n_data+i-1]);
    /*Segunda pendiente*/
    x_op=d[i-1]+h/2.0;
    y_op=d[n_data+i-1]+(h/2.0)*k1;
    k2=func(x_op,y_op);
    /*Tercera pendiente*/
    x_op=d[i-1]+h/2.0;
    y_op=d[n_data+i-1]+(h/2.0)*k2;
    k3=func(x_op, y_op);
    /*Cuarta Pendiente*/
    x_op=d[i-1]+h;
    y_op=d[n_data+i-1]+h*k3;
    k4=func(x_op, y_op);
    /*Calculo k*/
    k=(1.0/6.0)*(k1+(2.0*k2)+(2.0*k3)+k4);
    
    d[i]=d[i-1]+h;
    d[n_data+i]=d[n_data+i-1]+h*k;
    printf("%f\t%f\n",d[i], d[n_data+i]);
  }

  /* printf("%d\n",n_data);*/

  return 0;

}
