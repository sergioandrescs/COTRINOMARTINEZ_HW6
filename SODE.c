#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MSOLAR 1.989E30
#define UA 1.496E11
#define YEAR_SEC 31536000
#define G 6.67E-11
#define G_GRAV 39.486 //units of ua+3 msun-1 yr-1


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
  double k1_p1, k2_p1, k3_p1, k4_p1, ak1;
  double k1_p2, k2_p2, k3_p2, k4_p2, ak2;
  int n_data,i;
  double *d;  
  double x_op, y1_op, y2_op; 
  double new_G =(G*MSOLAR*pow(YEAR_SEC,2))/pow(UA,3);

  x_max=6.0;
  x_min=0.0;
  h=0.01;
  n_data=((x_max-x_min)/h);

  d=malloc(n_data*3*sizeof(double));

  d[0]=x_min;
  d[n_data]=1.0;
  d[2*n_data]=0;

  double fun1 (double x, double y1, double y2){
    return y2 ;
  }
  double fun2 (double x, double y1, double y2 ){
    return  -G_GRAV *mass[j]/ pow(r_ij,1.5) * (x[i] - x[j]);
  }


  for (i=1;i<n_data;i++){
    /*Primera pendiente*/   
    k1_p1=fun1(d[0],d[n_data+i-1],d[(n_data*2)+i-1]);
    k1_p2=fun2(d[0],d[n_data+i-1],d[(n_data*2)+i-1]);
    /*Segunda pendiente*/
    x_op=d[i-1]+h/2.0;
    y1_op=d[n_data+i-1]+(h/2.0)*k1_p1;
    y2_op=d[(2*n_data)+i-1]+(h/2.0)*k1_p2;
    k2_p1=fun1(x_op,y1_op,y2_op);
    k2_p2=fun2(x_op,y1_op,y2_op);
    /*Tercera pendiente*/
    x_op=d[i-1]+h/2.0;
    y1_op=d[n_data+i-1]+(h/2.0)*k2_p1;
    y2_op=d[(2*n_data)+i-1]+(h/2.0)*k2_p2;
    k3_p1=fun1(x_op,y1_op,y2_op);
    k3_p2=fun2(x_op,y1_op,y2_op);
    /*Cuarta Pendiente*/
    x_op=d[i-1]+h;
    y1_op=d[n_data+i-1]+(h)*k3_p1;
    y2_op=d[(2*n_data)+i-1]+(h)*k3_p2;
    k4_p1=fun1(x_op,y1_op,y2_op);
    k4_p2=fun2(x_op,y1_op,y2_op);
    /*Calculo k*/
    ak1=(1.0/6.0)*(k1_p1+(2.0*k2_p1)+(2.0*k3_p1)+k4_p1);
    ak2=(1.0/6.0)*(k1_p2+(2.0*k2_p2)+(2.0*k3_p2)+k4_p2);

    d[i]=d[i-1]+h;
    d[n_data+i]=d[n_data+i-1]+h*ak1;
    d[(2*n_data)+i]=d[(2*n_data)+i-1]+h*ak2;
    printf("%f\t%f\n",d[i], d[n_data+i]);
  }

  /* printf("%d\n",n_data);*/

  return 0;

}
