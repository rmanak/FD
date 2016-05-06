#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <bbhutil.h>
#include <init_f.h>
#include <residual_f.h>
#include <update_f.h>
/* Shapes: */
int Nx;
int shape[1];
int dim;
int level;
/* Coordinates: */
double *x;
/* Grid Functions: */
double *n_f;
double *nm1_f;
double *np1_f;
/* Parameters: */
double T0;
double T1;
double ht;
double myzero;
double xmax;
double xmin;
/* Coordinate Parameters: */
double x_max;
double x_min;
double hx;
double bbox[2];
int phys_bdy[2];
/* Time Evolution Parameters: */
int steps;
int output_freq;
double lambda;
double time;
void swap_levels( double **a, double **b){
	double *t;
	t = *a;
	*a = *b;
	*b = t;
}
void read_params(char *p_file,double *T1,double *lambda,int *Nx,double *x_max,double *myzero,double *T0,double *xmax,int *output_freq,int *level,double *x_min,double *xmin,int *steps)
{
get_param(p_file,"T1","double",1,T1);
get_param(p_file,"lambda","double",1,lambda);
get_param(p_file,"Nx","long",1,Nx);
get_param(p_file,"x_max","double",1,x_max);
get_param(p_file,"myzero","double",1,myzero);
get_param(p_file,"T0","double",1,T0);
get_param(p_file,"xmax","double",1,xmax);
get_param(p_file,"output_freq","long",1,output_freq);
get_param(p_file,"level","long",1,level);
get_param(p_file,"x_min","double",1,x_min);
get_param(p_file,"xmin","double",1,xmin);
get_param(p_file,"steps","long",1,steps);
}

int main(int argc, char **argv) {
char pfile[64];
strcpy(pfile,argv[1]);
/* Initialization of Coordinate: */
dim =1;
read_params(pfile,&T1,&lambda,&Nx,&x_max,&myzero,&T0,&xmax,&output_freq,&level,&x_min,&xmin,&steps);

Nx = Nx*(int)pow(2.0,(double)level)+1;
steps = steps*(int)pow(2.0,(double)level);
/* Allocating Memory to Grid Functions: */
x = vec_alloc(Nx);
n_f = vec_alloc(1*Nx);
nm1_f = vec_alloc(1*Nx);
np1_f = vec_alloc(1*Nx);
hx = (x_max-x_min)/(Nx-1);
dvumsh(x,Nx,x_min,x_max);
ht = lambda*hx*hx;
shape[0]=Nx;
bbox[0]=x_min;
bbox[1]=x_max;
time=0.0;
init_f_(x,&Nx,&T0,&T1,&xmax,&xmin,n_f);
gft_out_bbox("f.sdf",time,shape,dim,bbox,n_f);
rdvcpy(np1_f,n_f,Nx);
phys_bdy[0] = 1;
phys_bdy[1] = 1;
int i;
double tol,tres,norm_f;
for (i=0; i<steps; i++) {
   norm_f = l2norm(Nx,n_f);
   update_f_(n_f,np1_f,x,&Nx,&T0,&T1,&ht,&hx,&myzero,phys_bdy,np1_f);
   residual_f_(n_f,np1_f,x,&Nx,&T0,&T1,&ht,&hx,&myzero,phys_bdy,&tres);
   tol = tres/norm_f;
	time = time + ht;
	if ((i + 1) % (output_freq*(int)pow(2.0,(double)level))  == 0) {
	        gft_out_bbox("f.sdf",time,shape,dim,bbox,np1_f);
           printf("step: %d time : %f res: %1.14e\n",i+1,time,tol);	
        }
    swap_levels(&np1_f,&n_f);
  }
gft_close_all();
}
