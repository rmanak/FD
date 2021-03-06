
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <bbhutil.h>
#include <sdf_read.h>
#include <u_f.h>
#include <init_f.h>
#include <init_v.h>
#include <ire_f.h>
#include <res_f.h>
#include <res_f_t.h>
#include <u_f_t.h>
#include <init_f_t.h>
/* Shapes: */
int Nx;
int shape[1];
int dim;
int level;
/* Coordinates: */
double *x;
/* Grid Functions: */
double *n_f;
double *v;
double *n_f_t;
double *nm1_f;
double *nm1_f_t;
double *np1_f;
double *np1_f_t;
double myzero;
/* Parameters: */
double A;
double a;
double b;
double d;
double ht;
double r;
double x0;
double delx;
double idsignum;
double res_f,res_f_t;
double vcons;
double p;
int output_freq;
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

double tot_res;
void swap_levels(double **a, double **b) {
     double *t;
     t = *a;
     *a = *b;
     *b = t;
}

void copy_lev(double **a, double **b) {
	*b = *a;
}



void read_params(char *p_file,double *x_min,double *d,double *A,double *delx,double *r,int *steps,double *lambda,double *myzero,double *b,double *x_max,int *level,double *a,double *x0,int *Nx,double *idsignum, int *output_freq,double *vcons, double *p)
{
get_param(p_file,"x_min","double",1,x_min);
get_param(p_file,"d","double",1,d);
get_param(p_file,"A","double",1,A);
get_param(p_file,"delx","double",1,delx);
get_param(p_file,"r","double",1,r);
get_param(p_file,"steps","long",1,steps);
get_param(p_file,"lambda","double",1,lambda);
get_param(p_file,"myzero","double",1,myzero);
get_param(p_file,"b","double",1,b);
get_param(p_file,"x_max","double",1,x_max);
get_param(p_file,"level","long",1,level);
get_param(p_file,"a","double",1,a);
get_param(p_file,"x0","double",1,x0);
get_param(p_file,"Nx","long",1,Nx);
get_param(p_file,"output_freq","long",1,output_freq);
get_param(p_file,"idsignum","double",1,idsignum);
get_param(p_file,"vcons","double",1,vcons);
get_param(p_file,"p","double",1,p);
}

int main(int argc, char **argv) {
int j;
char pfile[64];
strcpy(pfile,argv[1]);
/* Initialization of Coordinate: */
dim =1;
read_params(pfile,&x_min,&d,&A,&delx,&r,&steps,&lambda,&myzero,&b,&x_max,&level,&a,&x0,&Nx,&idsignum,&output_freq,&vcons,&p);

Nx = Nx*(int)pow(2.0,(double)level)+1;
steps = steps*(int)pow(2.0,(double)level);
/* Allocating Memory to Grid Functions: */
x = vec_alloc(Nx);
n_f = vec_alloc(1*Nx);
v = vec_alloc(1*Nx);
n_f_t = vec_alloc(1*Nx);
nm1_f = vec_alloc(1*Nx);
nm1_f_t = vec_alloc(1*Nx);
np1_f = vec_alloc(1*Nx);
np1_f_t = vec_alloc(1*Nx);
hx = (x_max-x_min)/(Nx-1);
dvumsh(x,Nx,x_min,x_max);
ht = lambda*sqrt( 0.0+ hx*hx);
shape[0]=Nx;
bbox[0]=x_min;
bbox[1]=x_max;
time=0.0;
init_f_(x,&Nx,&A,&delx,&x0,n_f);
init_f_(x,&Nx,&A,&delx,&x0,np1_f);

init_v_(x,&Nx,&d,&p,&vcons,v);

init_f_t_(v,x,&Nx,&A,&delx,&idsignum,&x0,n_f_t);
init_f_t_(v,x,&Nx,&A,&delx,&idsignum,&x0,np1_f_t);


gft_out_bbox("f.sdf",time,shape,dim,bbox,np1_f);
gft_out_bbox("f_t.sdf",time,shape,dim,bbox,np1_f_t);
printf("step: %d time: %f iter: %d res: %1.14e ire: %1.14e\n",0,time,0,0.0,0.0);

int i;

for (i=0; i<steps; i++) {
  j=0;
  tot_res = 1.0;
  ire_f = 0;
  while (tot_res > 1.0e-9) {
      u_f_(n_f,n_f_t,np1_f,np1_f_t,v,x,&Nx,&ht,&hx,&myzero,phys_bdy,np1_f);
      u_f_t_(n_f,n_f_t,np1_f,np1_f_t,v,x,&Nx,&a,&b,&ht,&hx,&myzero,phys_bdy,np1_f_t);
      res_f_(n_f,n_f_t,np1_f,np1_f_t,v,x,&Nx,&ht,&hx,&myzero,phys_bdy,&res_f);
      res_f_t_(n_f,n_f_t,np1_f,np1_f_t,v,x,&Nx,&a,&b,&ht,&hx,&myzero,phys_bdy,&res_f_t);
		tot_res=res_f/l2norm(Nx,np1_f) + res_f_t/l2norm(Nx,np1_f_t);
		j++;
      }

   if (i>1) {
     ire_f_(n_f,nm1_f,np1_f,v,&Nx,&a,&b,&ht,&hx,&ire_f);
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
