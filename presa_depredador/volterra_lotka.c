#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define A 20.0
#define B 1.0
#define C 30.0
#define D 1.0

int main(){

  /*Definicion variables*/
  int i;
  float h;
  float t_min;
  float t_max;
  float n_puntos;
  int puntos;
  float k1x;
  float k1y;
  float t1;
  float x1;
  float y1;
  float k2x;
  float k2y;
  float t2;
  float x2;
  float y2;
  float k3x;
  float k3y;
  float t3;
  float x3;
  float y3;
  float k4x;
  float k4y;
  float k_medio_x;
  float k_medio_y;
  FILE *in;
  
  /*Inicializacion variables*/
  h=0.001; //Paso de iteracion
  t_min=0.0;  //Rango de tiempo
  t_max=1.0;
  n_puntos=(t_max-t_min)/h;  //Numero de puntos
  puntos=(int)n_puntos+1;
  printf("%d\n",puntos);
  float t[puntos];  //Variable independiente
  float x[puntos];  //Variables dependientes
  float y[puntos];
  in=fopen("datos.txt","w");

  /*Funcion x'(t)*/
  float x_prima(T,X,Y){
    float x_punto;
    x_punto= A*X-B*X*Y;
    return x_punto;
  }

  /*Funcion y'(t)*/
  float y_prima(T,X,Y){
    float y_punto;
    y_punto= -C*Y+D*X*Y;
    return y_punto;
  }
  
  /*Condiciones iniciales X(0) y Y(0) para  X'(0)=Y'(0)=0*/
  t[0]=0;  //Tiempo inicial
  y[0]=A/B; //Al evaluar X'(0)=0 se puede cancelar X(0) en la ecuacion
  x[0]=C/D; //Al evaluar Y'(0)=0 se puede cancelar Y(0) en la ecuacion

  do{
  /*Solucion del conjunto de ecuaciones diferenciales variando el valor inicial de X*/
  for (i=1;i<puntos;i++){
      k1x=x_prima(t[i-1],x[i-1],y[i-1]);
      k1y=y_prima(t[i-1],x[i-1],y[i-1]);

      /*Primera iteracion*/
      t1=t[i-1]+(h/2.0);      
      x1=x[i-1]+(h/2.0)*k1x;
      y1=y[i-1]+(h/2.0)*k1y;
      k2x=x_prima(t1,x1,y1);
      k2y=y_prima(t1,x1,y1);

      /*Segunda iteracion*/
      t2=t[i-1]+(h/2.0);      
      x2=x[i-1]+(h/2.0)*k2x;
      y2=y[i-1]+(h/2.0)*k2y;
      k3x=x_prima(t2,x2,y2);
      k3y=y_prima(t2,x2,y2);
      
      /*Tercera iteracion*/
      t3=t[i-1]+h;      
      x3=x[i-1]+h*k3x;
      y3=y[i-1]+h*k3y;
      k4x=x_prima(t3,x3,y3);
      k4y=y_prima(t3,x3,y3);

      /*Cuarta iteracion*/
      k_medio_x=(1.0/6.0)*(k1x+2.0*k2x+2.0*k3x+k4x);
      k_medio_y=(1.0/6.0)*(k1y+2.0*k2y+2.0*k3y+k4y);
    
      t[i]=t[i-1]+h;
      x[i]=x[i-1]+h*k_medio_x;
      y[i]=y[i-1]+h*k_medio_y;
      fprintf(in,"%f %f\n",x[i],y[i]);
  }
  printf("%f\n",x[0]);
  x[0]=x[0]-1;
  }while(x[0]>0);

  fclose(in);

return 0;

}
