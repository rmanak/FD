
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
int Ny;
int shape[2];
int dim;
int level;
/* Coordinates: */
double *x;
double *y;
/* Grid Functions: */
double *n_f;
double *n_f_t;
double *nm1_f;
double *nm1_f_t;
double *np1_f;
double *np1_f_t;
double res_f,res_f_t,ire_f;
double tot_res;
/* Parameters: */
double A;
double ht;
double xc;
double yc;
double delx;
double dely;
double idsigx;
double idsigy;
double myzero;
/* Coordinate Parameters: */
double x_max;
double x_min;
double hx;
double y_max;
double y_min;
double hy;
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

void set_cons(double *f, int N, double c) {
	   int i;
		   for(i=0;i<N;i++) f[i] = c;
}


void read_params(char *p_file,int *steps,double *idsigy,double *x_max,int *level,double *idsigx,double *delx,int *output_freq,double *lambda,double *dely,double *xc,double *yc,int *Ny,double *myzero,double *x_min,double *A,double *y_min,double *y_max,int *Nx)
{
get_param(p_file,"steps","long",1,steps);
get_param(p_file,"idsigy","double",1,idsigy);
get_param(p_file,"x_max","double",1,x_max);
get_param(p_file,"level","long",1,level);
get_param(p_file,"idsigx","double",1,idsigx);
get_param(p_file,"delx","double",1,delx);
get_param(p_file,"output_freq","long",1,output_freq);
get_param(p_file,"lambda","double",1,lambda);
get_param(p_file,"dely","double",1,dely);
get_param(p_file,"xc","double",1,xc);
get_param(p_file,"yc","double",1,yc);
get_param(p_file,"Ny","long",1,Ny);
get_param(p_file,"myzero","double",1,myzero);
get_param(p_file,"x_min","double",1,x_min);
get_param(p_file,"A","double",1,A);
get_param(p_file,"y_min","double",1,y_min);
get_param(p_file,"y_max","double",1,y_max);
get_param(p_file,"Nx","long",1,Nx);
}

int main(int argc, char **argv) {
char pfile[64];
int j;
strcpy(pfile,argv[1]);
/* Initialization of Coordinate: */
dim =2;
read_params(pfile,&steps,&idsigy,&x_max,&level,&idsigx,&delx,&output_freq,&lambda,&dely,&xc,&yc,&Ny,&myzero,&x_min,&A,&y_min,&y_max,&Nx);
Nx = Nx*(int)pow(2.0,(double)level)+1;
steps = steps*(int)pow(2.0,(double)level);
Ny = Ny*(int)pow(2.0,(double)level)+1;
/* Allocating Memory to Grid Functions: */
x = vec_alloc(Nx);
y = vec_alloc(Ny);
n_f = vec_alloc(1*Nx*Ny);
n_f_t = vec_alloc(1*Nx*Ny);
nm1_f = vec_alloc(1*Nx*Ny);
nm1_f_t = vec_alloc(1*Nx*Ny);
np1_f = vec_alloc(1*Nx*Ny);
np1_f_t = vec_alloc(1*Nx*Ny);
hx = (x_max-x_min)/(Nx-1);
dvumsh(x,Nx,x_min,x_max);
hy = (y_max-y_min)/(Ny-1);
dvumsh(y,Ny,y_min,y_max);
ht = lambda*sqrt( 0.0+ hx*hx+ hy*hy);
shape[0]=Nx;
shape[1]=Ny;
bbox[0]=x_min;
bbox[1]=x_max;
bbox[2]=y_min;
bbox[3]=y_max;
phys_bdy[0] = 1;
phys_bdy[1] = 1;
phys_bdy[2] = 1;
phys_bdy[3] = 1;
time=0.0;
init_f_(x,y,&Nx,&Ny,&A,&xc,&yc,&delx,&dely,n_f);
init_f_(x,y,&Nx,&Ny,&A,&xc,&yc,&delx,&dely,np1_f);
init_f_t_(x,y,&Nx,&Ny,&A,&xc,&yc,&delx,&dely,&idsigx,&idsigy,n_f_t);
init_f_t_(x,y,&Nx,&Ny,&A,&xc,&yc,&delx,&dely,&idsigx,&idsigy,np1_f_t);
gft_out_bbox("f.sdf",time,shape,dim,bbox,np1_f);
gft_out_bbox("f_t.sdf",time,shape,dim,bbox,np1_f_t);
	
for (int i=0; i<steps; i++) {
		  
	  j=0;
	  tot_res = 1.0;
	  ire_f = 0;
	  while (tot_res > 1.0e-9) {
          u_f_(n_f,x,y,n_f_t,np1_f,np1_f_t,&Nx,&Ny,&ht,&myzero,phys_bdy,np1_f);
          u_f_t_(n_f,x,y,n_f_t,np1_f,np1_f_t,&Nx,&Ny,&ht,&hx,&hy,&myzero,phys_bdy,np1_f_t);
          res_f_(n_f,x,y,n_f_t,np1_f,np1_f_t,&Nx,&Ny,&ht,&myzero,phys_bdy,&res_f);
          res_f_t_(n_f,x,y,n_f_t,np1_f,np1_f_t,&Nx,&Ny,&ht,&hx,&hy,&myzero,phys_bdy,&res_f_t);
			 tot_res = res_f/l2norm(Nx*Ny,np1_f) + res_f_t/l2norm(Nx*Ny,np1_f_t);
          j++;
	     } //End of while
		  time = time + ht;
        if (i>1) {
          ire_f_(np1_f,n_f,nm1_f,&Nx,&Ny,&ht,&hx,&hy,&ire_f);
		  }
		  if ((i + 1) % (output_freq*(int)pow(2.0,(double)level))  == 0) {
           gft_out_bbox("f.sdf",time,shape,dim,bbox,np1_f);
   	     gft_out_bbox("f_t.sdf",time,shape,dim,bbox,np1_f_t);
		     printf("step: %d time: %f iter: %d res: %1.14e ire: %1.14e\n",i+1,time,j,tot_res,ire_f);
		  }
		  rdvcpy(nm1_f,n_f,Nx*Ny);
		  swap_levels(&np1_f,&n_f);
		  swap_levels(&np1_f_t,&n_f_t);
		  

   }//End of for
gft_close_all();
}
