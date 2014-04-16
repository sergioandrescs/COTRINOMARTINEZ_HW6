#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define G_GRAV 39.486 //units of ua+3 msun-1 yr-1
#define DELTA 0.01

/*Codigo que genera las condiciones iniciales de posicion y velocidad para una galaxia.

Para ejecutar el codigo correctamente, es necesario ingresar la informacion en el siguiente orden al momento de ejecutar el codigo:
Masa central, Radio de la galaxia, Numero de elemento, posicion inicial de la galaxia (x0, y0, z0), velocidad inicial de la galaxia (vx0, vy0, vz0)*/

FLOAT *get_memory(int n_points);
void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT x0, FLOAT y0, FLOAT z0, int n_points, FLOAT radius);
void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT vx0, FLOAT vy0, FLOAT vz0, int n_points, FLOAT vel, FLOAT radius);
void print_initial(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT *id);
void make_id(FLOAT *id, int n_points);

int main(int argc, char **argv){

  /*i for iteration*/
  int i;

  /*ID of all particles*/
  FLOAT *id;

  /*Initial Position*/
  FLOAT x0 = atof(argv[4]);
  FLOAT y0 = atof(argv[5]);
  FLOAT z0 = atof(argv[6]);

  /*Initial Velocity*/
  FLOAT vx0 = atof(argv[7]);
  FLOAT vy0 = atof(argv[8]);
  FLOAT vz0 = atof(argv[9]);

  /*positions of all particles*/
  FLOAT *x;
  FLOAT *y;
  FLOAT *z;
  
  /*velocities of all particles*/
  FLOAT *v_x;
  FLOAT *v_y;
  FLOAT *v_z;

  int n_points = atoi(argv[3])+1;
  FLOAT radius = atof(argv[2]);
  FLOAT unit_mass = atof(argv[1]);
  FLOAT vel_initial = 1; 

  /*memory allocation*/
  x = get_memory(n_points);
  y = get_memory(n_points);
  z = get_memory(n_points);
  v_x = get_memory(n_points);
  v_y = get_memory(n_points);
  v_z = get_memory(n_points);
  id = get_memory(n_points);

  make_id(id,n_points);
  initialize_pos(x,y,z,x0,y0,z0,n_points,radius);
  initialize_vel(v_x,v_y,v_z,vx0,vy0,vz0,n_points,vel_initial,radius);
  print_initial(x,y,z,v_x,v_y,v_z,n_points,id);

}

void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT x0, FLOAT y0, FLOAT z0, int n_points, FLOAT radius){
  int i; 
  FLOAT delta_theta;
  delta_theta = 2.0*PI/n_points;

  x[0]=x0;
  y[0]=y0;
  z[0]=z0;
  
  for(i=0;i<n_points-1;i++){
    x[i+1] = cos(delta_theta * i) * radius + x[0];
    y[i+1] = sin(delta_theta * i) * radius + y[0];
    z[i+1] = 0.0 + z[0];
  }
}

void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT vx0, FLOAT vy0, FLOAT vz0, int n_points, FLOAT vel, FLOAT radius){
  int i; 
  FLOAT delta_theta;
  delta_theta = 2.0*PI/n_points;

  vx[0] = vx0;
  vy[0] = vy0;
  vz[0] = vz0;
  
  for(i=0;i<n_points-1;i++){
    vx[i+1] = -sin(delta_theta * i) * vel + vx[0];
    vy[i+1] = cos(delta_theta * i) * vel + vy[0];
    vz[i+1] = 0.0 + vz[0];
  }  
}

void make_id(FLOAT *id, int n_points){
  int i;

  for (i=0;i<n_points;i++){
    if(i==0){
      id[i] = -1.0*(i+1);
    }
    else{
      id[i] = (i+1);
    }
  }
}


FLOAT * get_memory(int n_points){
  FLOAT * x; 
  if(!(x = malloc(sizeof(FLOAT) * n_points))){
    printf("problem with memory allocation");
    exit(1);
  }
  return x;
}

void print_initial(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT *id){
  int i;
  FILE *in;
  in = fopen("inicialdata.dat","w");
  for(i=0;i<n_points;i++){
    fprintf(in,"%f\t %f\t %f\t %f\t %f\t %f\t %f  \n", 
	   id[i], x[i], y[i], z[i], vx[i], vy[i], vz[i]);
  }
  fclose(in);
}

