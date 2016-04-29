
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <bbhutil.h>
#include <sdf_read.h>
#include <u_f.h>
#include <init_f.h>
#include <ire_f.h>
#include <res_f.h>
#include <res_f_t.h>
#include <u_f_t.h>
#include <init_f_t.h>
/* Shapes: */
int Nx;
int Nz;
int shape[2];
int dim;
int level;
/* Coordinates: */
double *x;
double *z;
/* Grid Functions: */
double *n_f_t;
double *np1_f_t;
double *n_f;
double *n_rf;
double *nm1_rf;
double *nm2_rf;
double *np1_f;
double *np1_rf;
double *np2_rf;
/* Parameters: */
double A;
double ht;
double xc;
double zc;
double delx;
double delz;
double idsigx;
double idsigz;
double myzero;
/* Coordinate Parameters: */
double x_max;
double x_min;
double hx;
double z_max;
double z_min;
double hz;
double bbox[4];
int phys_bdy[4];
/* Time Evolution Parameters: */
int steps;
int output_freq;
double lambda;
double time;
void swap_levels(double **a, double **b) {
	     double *t;
		       t = *a;
				      *a = *b;
						     *b = t;
}

void read_params(char *p_file,double *z_min,double *zc,double *delx,double *A,double *lambda,double *x_min,double *delz,int *Nz,double *x_max,double *z_max,double *idsigx,double *xc,double *myzero,int *steps,int *Nx,int *output_freq,double *idsigz,int *level)
{
get_param(p_file,"z_min","double",1,z_min);
get_param(p_file,"zc","double",1,zc);
get_param(p_file,"delx","double",1,delx);
get_param(p_file,"A","double",1,A);
get_param(p_file,"lambda","double",1,lambda);
get_param(p_file,"x_min","double",1,x_min);
get_param(p_file,"delz","double",1,delz);
get_param(p_file,"Nz","long",1,Nz);
get_param(p_file,"x_max","double",1,x_max);
get_param(p_file,"z_max","double",1,z_max);
get_param(p_file,"idsigx","double",1,idsigx);
get_param(p_file,"xc","double",1,xc);
get_param(p_file,"myzero","double",1,myzero);
get_param(p_file,"steps","long",1,steps);
get_param(p_file,"Nx","long",1,Nx);
get_param(p_file,"output_freq","long",1,output_freq);
get_param(p_file,"idsigz","double",1,idsigz);
get_param(p_file,"level","long",1,level);
}

int main(int argc, char **argv) {
char pfile[64];
int j;
double tot_res;
double res_f,res_f_t;
double ire_f;
double norm_f,norm_f_t;
strcpy(pfile,argv[1]);
/* Initialization of Coordinate: */
dim =2;
read_params(pfile,&z_min,&zc,&delx,&A,&lambda,&x_min,&delz,&Nz,&x_max,&z_max,&idsigx,&xc,&myzero,&steps,&Nx,&output_freq,&idsigz,&level);

Nx = Nx*(int)pow(2.0,(double)level)+1;
steps = steps*(int)pow(2.0,(double)level);
Nz = Nz*(int)pow(2.0,(double)level)+1;
/* Allocating Memory to Grid Functions: */
x = vec_alloc(Nx);
z = vec_alloc(Nz);
n_f = vec_alloc(1*Nx*Nz);
n_f_t = vec_alloc(1*Nx*Nz);
np1_f_t = vec_alloc(1*Nx*Nz);
nm1_rf = vec_alloc(1*Nx*Nz);
nm2_rf = vec_alloc(1*Nx*Nz);
n_rf = vec_alloc(1*Nx*Nz);
np1_rf = vec_alloc(1*Nx*Nz);
np1_f = vec_alloc(1*Nx*Nz);
np2_rf = vec_alloc(1*Nx*Nz);
hx = (x_max-x_min)/(Nx-1);
dvumsh(x,Nx,x_min,x_max);
hz = (z_max-z_min)/(Nz-1);
dvumsh(z,Nz,z_min,z_max);
ht = lambda*sqrt( 0.0+ hx*hx+ hz*hz);
shape[0]=Nx;
shape[1]=Nz;
bbox[0]=x_min;
bbox[1]=x_max;
bbox[2]=z_min;
bbox[3]=z_max;
time=0.0;
init_f_(z,x,&Nx,&Nz,&A,&xc,&zc,&delx,&delz,n_f);
init_f_(z,x,&Nx,&Nz,&A,&xc,&zc,&delx,&delz,np1_f);
init_f_t_(z,x,&Nx,&Nz,&A,&xc,&zc,&delx,&delz,&idsigx,&idsigz,n_f_t);
init_f_t_(z,x,&Nx,&Nz,&A,&xc,&zc,&delx,&delz,&idsigx,&idsigz,np1_f_t);

gft_out_bbox("f.sdf",time,shape,dim,bbox,np1_f);
gft_out_bbox("f_t.sdf",time,shape,dim,bbox,np1_f_t);


for (int i=0; i<steps; i++) {
    tot_res = 1.0;
	 j = 0;
	 ire_f = 0;
	 while ( tot_res > 1.0e-9) {
     u_f_(n_f,x,n_f_t,np1_f,np1_f_t,&Nx,&Nz,&ht,&hx,&hz,phys_bdy,np1_f);
     u_f_t_(n_f,x,n_f_t,np1_f,np1_f_t,&Nx,&Nz,&ht,&hx,&hz,phys_bdy,np1_f_t);
     res_f_(n_f,x,n_f_t,np1_f,np1_f_t,&Nx,&Nz,&ht,&hx,&hz,phys_bdy,&res_f);
     res_f_t_(n_f,x,n_f_t,np1_f,np1_f_t,&Nx,&Nz,&ht,&hx,&hz,phys_bdy,&res_f_t);
     norm_f = 1.0/l2norm(Nx*Nz,np1_f);
     norm_f_t = 1.0/l2norm(Nx*Nz,np1_f_t);
	  j++;
	  tot_res = res_f*norm_f + res_f_t*norm_f_t;
	  if (j>300) {
		  printf("iteration did not converge\n");
		  printf("res = %f\n",res_f+res_f_t);
		  exit(1);
	  }
	 }

	  if (i>5){ 
      ire_f_(n_rf,nm2_rf,np1_rf,nm1_rf,x,np2_rf,&Nx,&Nz,&ht,&hx,&hz,&ire_f);
	  }

     time = time + ht;

	  if ((i + 1) % (output_freq*(int)pow(2.0,(double)level))  == 0) {
        gft_out_bbox("f.sdf",time,shape,dim,bbox,np1_f);
        gft_out_bbox("f_t.sdf",time,shape,dim,bbox,np1_f_t);
		  printf("step: %d time: %f iter: %d res: %1.14e ire: %1.14e\n",i+1,time,j,tot_res,ire_f);
	  }

	  rdvcpy(nm2_rf,nm1_rf,Nx*Nz);
	  rdvcpy(nm1_rf,n_rf,Nx*Nz);
	  rdvcpy(n_rf,np1_rf,Nx*Nz);
	  rdvcpy(np1_rf,np2_rf,Nx*Nz);
	  rdvcpy(np2_rf,np1_f,Nx*Nz);

	  swap_levels(&np1_f,&n_f);
	  swap_levels(&np1_f_t,&n_f_t);

}
}
//init_f_(z,x,&Nx,&Nz,&A,&xc,&zc,&delx,&delz,res);
//init_f_t_(z,x,&Nx,&Nz,&A,&xc,&zc,&delx,&delz,&idsigx,&idsigz,res);

