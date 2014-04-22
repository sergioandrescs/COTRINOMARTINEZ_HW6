#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define G_GRAV 39.486 //units of ua+3 msun-1 yr-1
#define DELTA 0.01

/*El siguiente codigo modela la colision de dos galaxias consistentes en una masa central y una cantidad de particulas que giran a su alrededor.

El codigo recibe, en ese orden: nombre del archivo con datos iniciales, el tiempo de evolucion T (en miles de millones de años), y el numero de momentos diferentes equiespaciados M.  */

FLOAT *get_memory(int n_points);
void print_status(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points, FLOAT energy);
void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points, FLOAT *id);
void runge_kutta4(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points, FLOAT *mass, FLOAT *id);
FLOAT get_energy (FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, FLOAT *mass);
void get_initial (FLOAT *id, FLOAT *x, FLOAT *y, FLOAT*z, FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, char *fileic);

int main(int argc, char **argv){
  /*Name initial data file*/
  char *fileic = (argv[1]);

  /*energy*/
  FLOAT energy;

  /*Momentos de momentos equiespaciados*/
  int momentos = atoi(argv[3]);

  /*i for iteration*/
  int i;

  /*ID of all particles*/
  FLOAT *id;

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
  int n_steps = (int)((1E9*atoi(argv[2]))/delta_t);
  int n_points;
  FLOAT radius = 100.0;
  FLOAT unit_mass = 5.0; 
  FLOAT vel_initial = 1.1;

  
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
  id = get_memory(n_points);

  get_acceleration(a_x,a_y,a_z,x,y,z,mass,n_points,id);
  energy = get_energy(x,y,z,v_x,v_y,v_z,n_points,mass);
  for (i=0;i<n_steps;i++){
    if(i==0){
      print_status(x,y,z,v_x,v_y,v_z, a_x, a_y, a_z, n_points,energy);
    }
    else{
      runge_kutta4(x,y,z,v_x,v_y,v_z,a_x,a_y,a_z,n_points,mass,id);
      energy = get_energy(x,y,z,v_x,v_y,v_z,n_points,mass);
      print_status(x,y,z,v_x,v_y,v_z, a_x, a_y, a_z, n_points,energy);
    }
  }
}

void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points, FLOAT *id){
  int i,j;
  FLOAT r_ij;
  int count;
  for(j=0;j<n_points;j++){
    if(id[j]<0){
      for(i=0;i<n_points;i++){
	ax[i]=0.0;
	ay[i]=0.0;
	az[i]=0.0;
	
	r_ij = (pow((x[i] - x[j]),2.0) +
		pow((y[i] - y[j]),2.0) +
		pow((z[i] - z[j]),2.0));
	r_ij = sqrt(r_ij);
	ax[i] += -G_GRAV *mass[0]/(0.0001 +  pow(r_ij,3)) * (x[i] - x[0]);
	ay[i] += -G_GRAV *mass[0]/(0.0001 +  pow(r_ij,3)) * (y[i] - y[0]);
	az[i] += -G_GRAV *mass[0]/(0.0001 +  pow(r_ij,3)) * (z[i] - z[0]);
      }
  } 
  }   
}

void get_initial (FLOAT *id, FLOAT *x, FLOAT *y, FLOAT*z, FLOAT *vx, FLOAT *vy, FLOAT *vz, int n_points, char *fileic){
  FILE *in;
  int i ,c;

  in = fopen(fileic, "r");
  do{
  c = fgetc(in);
  if(c=='\n'){
    n_points += 1;
  }
  
  }while(c!=EOF);

  if(!in){
    printf("Hay un problema con el archivo.\n");
    exit(1);
  }

  rewind(in);

   for(i=0;i<n_points*7;i++){ 
     fscanf(in, "%f %f %f %f %f %f %f \n", &(id[i]),&(x[i]), &(y[i]), &(z[i]), &(vx[i]), &(vy[i]), &(vz[i]) );
  }

 fclose(in);
}

FLOAT * get_memory(int n_points){
  FLOAT * x; 
  if(!(x = malloc(sizeof(FLOAT) * n_points))){
    printf("problem with memory allocation");
    exit(1);
  }
  return x;
}

void print_status(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points, FLOAT energy){
  int i;
  for(i=0;i<n_points;i++){
    printf("%f %f %f %f %f %f %f %f %f %f \n", 
	   x[i], y[i], z[i], energy, vx[i], vy[i], vz[i], ax[i], ay[i], az[i]);
  }
}
void runge_kutta4(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points, FLOAT *mass, FLOAT *id){
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
    get_acceleration( ax, ay, az, xts, yts, zts, mass, n_points, id);
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
    get_acceleration( ax, ay, az, xts, yts, zts, mass, n_points, id);
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
    get_acceleration( ax, ay, az, xts, yts, zts, mass, n_points, id);
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

