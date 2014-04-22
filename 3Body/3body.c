#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define G_GRAV 39.486 //units of ua+3 msun-1 yr-1
#define DELTA 0.01

/*El siguiente codigo modela el movimiento de 3 cuerpos con diferentes masas debido a sus campos gravitatorios. 
Como condiciones inicuales, se considera que la energia total del sistema (K+U) es igual a cero, que las pociciones iniciales son puntos equiespaciados sobre un circulo de radio dado,  y que la velocidad inicial de los cuerpos 
es perpendicular al plano en el que se encuentran dichos cuerpos.
El codigo imprime de vuelta las posiciones de los cuerpos en los 3 ejes, la energia total del sistema, la velocidad de los cuerpos y su aceleracion en los 3 ejes.*/

/*Inicializacion de funciones*/

FLOAT *get_memory(int n_points);
void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z, int n_points, FLOAT radius);
void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT vel, FLOAT radius);
void initialize_mass(FLOAT *mass, int n_points, FLOAT unit_mass);
void print_status(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points, FLOAT energy);
void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points);
void runge_kutta4(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points, FLOAT *mass);
FLOAT get_energy (FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT *mass);

int main(int argc, char **argv){
  /*energy*/
  FLOAT energy;

  /*i for iteration*/
  int i;

  /*positions of all particles*/
  FLOAT *x;
  FLOAT *y;
  FLOAT *z;
  
  /*velocities of all particles*/
  FLOAT *v_x;
  FLOAT *v_y;
  FLOAT *v_z;

  /*accelerations of all particles*/
  FLOAT *a_x;
  FLOAT *a_y;
  FLOAT *a_z;

  /*masses*/
  FLOAT *mass;

  /*timestep variables*/
  FLOAT delta_t= DELTA;
  int n_steps = (int)(100.0/delta_t);
  int n_points = 3;
  FLOAT radius = 100.0;
  FLOAT unit_mass = 1.0; 
  FLOAT vel_initial = sqrt((11.0/3.0) * G_GRAV * unit_mass / (sqrt(3.0)*radius));

  
  /*memory allocation*/
  x = get_memory(n_points);
  y = get_memory(n_points);
  z = get_memory(n_points);
  v_x = get_memory(n_points);
  v_y = get_memory(n_points);
  v_z = get_memory(n_points);
  a_x = get_memory(n_points);
  a_y = get_memory(n_points);
  a_z = get_memory(n_points);
  mass = get_memory(n_points);

  /*Se inicializan la posicion, velocidad, masas, aceleracion y energia
    Ademas, se crea un loop que recalcule la posicion de las masas usando el metodo Runge Kutta de orden 4*/

  initialize_pos(x,y,z, n_points, radius);
  initialize_vel(v_x,v_y,v_z, n_points, vel_initial, radius);
  initialize_mass(mass, n_points, unit_mass);
  get_acceleration(a_x,a_y,a_z,x,y,z,mass,n_points);
  energy = get_energy(x,y,z,v_x,v_y,v_z,n_points,mass);
  for (i=0;i<n_steps;i++){
    if(i==0){
      print_status(x,y,z,v_x,v_y,v_z, a_x, a_y, a_z, n_points,energy);
    }
    else{
      runge_kutta4(x,y,z,v_x,v_y,v_z,a_x,a_y,a_z,n_points,mass);
      energy = get_energy(x,y,z,v_x,v_y,v_z,n_points,mass);
      print_status(x,y,z,v_x,v_y,v_z, a_x, a_y, a_z, n_points,energy);
    }
  }
}

/*La aceleracion de cada masa se calcula en cada eje, usando la ley de la gravitacion universal y la segunda ley de Newton*/

void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points){
  int i,j;
  FLOAT r_ij;
  for(i=0;i<n_points;i++){
    ax[i]=0.0;
    ay[i]=0.0;
    az[i]=0.0;
    
    for(j=0;j<n_points;j++){
      if(j!=i){
	r_ij = (pow((x[i] - x[j]),2.0) +
		pow((y[i] - y[j]),2.0) +
		pow((z[i] - z[j]),2.0));
	r_ij = sqrt(r_ij);
	ax[i] += -G_GRAV *mass[j]/ pow(r_ij,1.5) * (x[i] - x[j]);
	ay[i] += -G_GRAV *mass[j]/ pow(r_ij,1.5) * (y[i] - y[j]);
	az[i] += -G_GRAV *mass[j] / pow(r_ij,1.5) * (z[i] - z[j]);
      }
    }    
  }  
}

/*Las pociciones iniciales son puntos equiespaciados sobre un circulo de un radio dado*/

void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z, int n_points, FLOAT radius){
  int i; 
  FLOAT delta_theta;
  delta_theta = 2.0*PI/n_points;
  
  for(i=0;i<n_points;i++){
    x[i] = cos(delta_theta * i) * radius;
    y[i] = sin(delta_theta * i) * radius;
    z[i] = 0.0;
  }
}

/*Las velocidades iniciales son perpendiculares al plano sobre el cual estan ubicadas las masas.*/

void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT vel, FLOAT radius){
  int i; 
  
  for(i=0;i<n_points;i++){
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = vel;
  }  

}

/*Las masas de los cuerpos estan medidas en masas solares. Cada masa es una unidad mayor que la anterior.*/

void initialize_mass(FLOAT *mass, int n_points, FLOAT unit_mass){
  int i;
  for (i=0;i<n_points;i++){
    mass[i] = (i+1) * unit_mass;
  }
}

/*Se realiza la reserva de memoria para los diferentes datos.*/

FLOAT * get_memory(int n_points){
  FLOAT * x; 
  if(!(x = malloc(sizeof(FLOAT) * n_points))){
    printf("problem with memory allocation");
    exit(1);
  }
  return x;
}

/*El codigo imprime de vuelta las posiciones de los cuerpos en los 3 ejes, la energia total del sistema, la velocidad de los cuerpos y su aceleracion en los 3 ejes*/

void print_status(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points, FLOAT energy){
  int i;
  for(i=0;i<n_points;i++){
    printf("%f %f %f %f %f %f %f %f %f %f \n", 
	   x[i], y[i], z[i], energy, vx[i], vy[i], vz[i], ax[i], ay[i], az[i]);
  }
}

/*Las posiciones de las masas se actualizan usando un metodo Runge Kutta de cuarto orden.*/

void runge_kutta4(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points, FLOAT *mass){
  int i; 
 
  /*slope for Runge Kutta Method*/
  FLOAT *kvx;
  FLOAT *kvy;
  FLOAT *kvz;
  FLOAT *kax;
  FLOAT *kay;
  FLOAT *kaz;

  /*test positions for Runge Kutta Method*/
  FLOAT *xts;
  FLOAT *yts;
  FLOAT *zts;

  /*timestep variables*/
  FLOAT delta_t= DELTA;

  /*memory allocation*/
  kvx = get_memory(4*n_points);
  kvy = get_memory(4*n_points);
  kvz = get_memory(4*n_points);
  kax = get_memory(4*n_points);
  kay = get_memory(4*n_points);
  kaz = get_memory(4*n_points);
  xts = get_memory(n_points);
  yts = get_memory(n_points);
  zts = get_memory(n_points);

  for(i=0;i<n_points;i++){
      /*Primera pendiente*/   
      kvx[(4*i)]=vx[i];
      kvy[(4*i)]=vy[i];
      kvz[(4*i)]=vz[i];
      kax[(4*i)]=ax[i];
      kay[(4*i)]=ay[i];
      kaz[(4*i)]=az[i];
    }
    for (i=0;i<n_points;i++){
      /*Segunda pendiente*/
      xts[i]=x[i]+(delta_t/2.0)*kvx[(4*i)];
      yts[i]=y[i]+(delta_t/2.0)*kvy[(4*i)];
      zts[i]=z[i]+(delta_t/2.0)*kvz[(4*i)];
      kvx[1+(4*i)]=vx[i]+(delta_t/2.0)*kax[(4*i)];
      kvy[1+(4*i)]=vy[i]+(delta_t/2.0)*kay[(4*i)];
      kvz[1+(4*i)]=vz[i]+(delta_t/2.0)*kaz[(4*i)];
    }
    get_acceleration( ax, ay, az, xts, yts, zts, mass, n_points);
    for(i=0;i<n_points;i++){
      kax[1+(4*i)]=ax[i];
      kay[1+(4*i)]=ay[i];
      kaz[1+(4*i)]=az[i];
    }

    for(i=0;i<n_points;i++){
      /*Tercera pendiente*/
      xts[i]=x[i]+(delta_t/2.0)*kvx[1+(4*i)];
      yts[i]=y[i]+(delta_t/2.0)*kvy[1+(4*i)];
      zts[i]=z[i]+(delta_t/2.0)*kvz[1+(4*i)];
      kvx[2+(4*i)]=vx[i]+(delta_t/2.0)*kax[1+(4*i)];
      kvy[2+(4*i)]=vy[i]+(delta_t/2.0)*kay[1+(4*i)];
      kvz[2+(4*i)]=vz[i]+(delta_t/2.0)*kaz[1+(4*i)];
    }
    get_acceleration( ax, ay, az, xts, yts, zts, mass, n_points);
    for(i=0;i<n_points;i++){
      kax[2+(4*i)]=ax[i];
      kay[2+(4*i)]=ay[i];
      kaz[2+(4*i)]=az[i];
    }

    for(i=0;i<n_points;i++){
      /*Cuarta Pendiente*/
      xts[i]=x[i]+(delta_t)*kvx[2+(4*i)];
      yts[i]=y[i]+(delta_t)*kvy[2+(4*i)];
      zts[i]=z[i]+(delta_t)*kvz[2+(4*i)];
      kvx[3+(4*i)]=vx[i]+(delta_t)*kax[2+(4*i)];
      kvy[3+(4*i)]=vy[i]+(delta_t)*kay[2+(4*i)];
      kvz[3+(4*i)]=vz[i]+(delta_t)*kaz[2+(4*i)];
    }
    get_acceleration( ax, ay, az, xts, yts, zts, mass, n_points);
    for(i=0;i<n_points;i++){
      kax[3+(4*i)]=ax[i];
      kay[3+(4*i)]=ay[i];
      kaz[3+(4*i)]=az[i];
    }

    for(i=0;i<n_points;i++){
      /*Calculo k*/
      vx[(i)]=(1.0/6.0)*(kvx[(4*i)]+(2.0*kvx[1+(4*i)])+(2.0*kvx[2+(4*i)])+kvx[3+(4*i)]);
      vy[(i)]=(1.0/6.0)*(kvy[(4*i)]+(2.0*kvy[1+(4*i)])+(2.0*kvy[2+(4*i)])+kvy[3+(4*i)]);
      vz[(i)]=(1.0/6.0)*(kvz[(4*i)]+(2.0*kvz[1+(4*i)])+(2.0*kvz[2+(4*i)])+kvz[3+(4*i)]);
      ax[(i)]=(1.0/6.0)*(kax[(4*i)]+(2.0*kax[1+(4*i)])+(2.0*kax[2+(4*i)])+kax[3+(4*i)]);
      ay[(i)]=(1.0/6.0)*(kay[(4*i)]+(2.0*kay[1+(4*i)])+(2.0*kay[2+(4*i)])+kay[3+(4*i)]);
      az[(i)]=(1.0/6.0)*(kaz[(4*i)]+(2.0*kaz[1+(4*i)])+(2.0*kaz[2+(4*i)])+kaz[3+(4*i)]);
    }
    for(i=0;i<n_points;i++){
      x[i]+=delta_t*vx[(i)];
      y[i]+=delta_t*vy[(i)];
      z[i]+=delta_t*vz[(i)];   
    }  
}
/*Se calcula la energia del sistema.*/

FLOAT get_energy (FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT *mass){
  int i,j;
  FLOAT r_ij;
  FLOAT energy = 0.0;
  for(i=0;i<n_points;i++){
    energy +=(0.5)*mass[i]*(pow(vx[i],2.0)+pow(vy[i],2.0)+pow(vz[i],2.0));
  }
  for(i=0;i<n_points-1;i++){    
    for(j=0;j<n_points;j++){
      if(j>i){
	r_ij = (pow((x[i] - x[j]),2.0) +
		pow((y[i] - y[j]),2.0) +
		pow((z[i] - z[j]),2.0));
	r_ij = sqrt(r_ij);
	energy += -G_GRAV *mass[j]*mass[i]/r_ij;
      }
    }    
  }  
  return energy;
}

