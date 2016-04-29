#include <stdlib.h>
#include <sdf_read.h>
#include <stdio.h>
#include <math.h>
#include <u_f.h>
#include <init_f.h>
#include <u_f_t.h>
#include <init_f_t.h>
#include <res_f.h>
#include <res_f_t.h>
#include <ire_f.h>
/* Shapes: */
int Nx;
int shape[1];
int dim;
int level;
/* Coordinates: */
double *x;
/* Grid Functions: */
double *n_f;
double *n_f_t;
double *np1_f;
double *np1_f_t;
double *nm1_f;
/* Parameters: */
double A;
double tot_res;
double ht;
double x0;
double delx;
double init_v;
double idsignum;
double res_f,res_f_t;
double ire_f;
/* Coordinate Parameters: */
double x_max;
double x_min;
double hx;
double bbox[2];
int phys_bdy[4];
/* Time Evolution Parameters: */
int steps;
double lambda;
double time;
int output_freq;
void swap_levels(double **a, double **b) {
     double *t;
     t = *a;
     *a = *b;
     *b = t;
}

void set_cons(double *f, int N, double c) {
	   int i;
		   for(i=0;i<N;i++) f[i] = c;
}

void copy_lev(double **a, double **b) {
	   *b = *a;
}


void read_params(char *p_file,double *A,double *x_min,double *lambda,double *idsignum,double *init_v,double *delx,int *level,double *x_max,int *Nx,int *steps,double *x0,int *output_freq)
{
get_param(p_file,"A","double",1,A);
get_param(p_file,"x_min","double",1,x_min);
get_param(p_file,"lambda","double",1,lambda);
get_param(p_file,"idsignum","double",1,idsignum);
get_param(p_file,"init_v","double",1,init_v);
get_param(p_file,"delx","double",1,delx);
get_param(p_file,"level","long",1,level);
get_param(p_file,"x_max","double",1,x_max);
get_param(p_file,"Nx","long",1,Nx);
get_param(p_file,"steps","long",1,steps);
get_param(p_file,"x0","double",1,x0);
get_param(p_file,"output_freq","long",1,output_freq);
}

int main(int argc, char **argv) {
char pfile[64];
int j;
strcpy(pfile,argv[1]);
/* Initialization of Coordinate: */
dim =1;
read_params(pfile,&A,&x_min,&lambda,&idsignum,&init_v,&delx,&level,&x_max,&Nx,&steps,&x0,&output_freq);

Nx = Nx*(int)pow(2.0,(double)level)+1;
steps = steps*(int)pow(2.0,(double)level);
/* Allocating Memory to Grid Functions: */
x = vec_alloc(Nx);
n_f = vec_alloc(1*Nx);
n_f_t = vec_alloc(1*Nx);
np1_f = vec_alloc(1*Nx);
np1_f_t = vec_alloc(1*Nx);
nm1_f = vec_alloc(Nx);
hx = (x_max-x_min)/(Nx-1);
dvumsh(x,Nx,x_min,x_max);
ht = lambda*sqrt( 0.0+ hx*hx);
shape[0]=Nx;
bbox[0]=x_min;
bbox[1]=x_max;
time=0.0;
init_f_(x,&Nx,&A,&delx,&x0,n_f);
init_f_(x,&Nx,&A,&delx,&x0,np1_f);
init_f_t_(x,&Nx,&A,&delx,&idsignum,&x0,n_f_t);
init_f_t_(x,&Nx,&A,&delx,&idsignum,&x0,np1_f_t);

phys_bdy[0] = 1;
phys_bdy[1] = 1;

gft_out_bbox("f.sdf",time,shape,dim,bbox,np1_f);
gft_out_bbox("f_t.sdf",time,shape,dim,bbox,np1_f_t);
printf("step: %d time: %f iter: %d res: %1.14e ire: %1.14e\n",0,time,0,0.0,0.0);


int i;
for (i=0; i<steps; i++) {
  j=0;
  tot_res = 1.0;
  ire_f = 0;
  while (tot_res > 1.0e-9) {
      u_f_(n_f,n_f_t,np1_f,np1_f_t,&Nx,&ht,phys_bdy,np1_f);
      u_f_t_(n_f,n_f_t,np1_f,np1_f_t,&Nx,&ht,&hx,phys_bdy,np1_f_t);
      res_f_(n_f,n_f_t,np1_f,np1_f_t,&Nx,&ht,phys_bdy,&res_f);
      res_f_t_(n_f,n_f_t,np1_f,np1_f_t,&Nx,&ht,&hx,phys_bdy,&res_f_t);
		tot_res=res_f/l2norm(Nx,np1_f) + res_f_t/l2norm(Nx,np1_f_t);
		j++;
      }
   if (i>1) {
     ire_f_(n_f,nm1_f,np1_f,&Nx,&ht,&hx,&ire_f);
	}
   time = time + ht;
	if ((i + 1) % (output_freq*(int)pow(2.0,(double)level))  == 0) {
	  gft_out_bbox("f.sdf",time,shape,dim,bbox,np1_f);
	  gft_out_bbox("f_t.sdf",time,shape,dim,bbox,np1_f_t);
	  printf("step: %d time: %f iter: %d res: %1.14e ire: %1.14e\n",i+1,time,j,tot_res,ire_f);
	}
	rdvcpy(nm1_f,n_f,Nx);
	swap_levels(&np1_f,&n_f);
	swap_levels(&np1_f_t,&n_f_t);
}
gft_close_all();
}
