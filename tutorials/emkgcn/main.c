
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <bbhutil.h>
#include <sdf_read.h>
#include <init_pi.h>
#include <init_pp.h>
#include <ire_pi.h>
#include <ire_pp.h>
#include <res_pi.h>
#include <res_pp.h>
#include <u_pi.h>
#include <u_pp.h>
#include <calc_trr.h>
#include <calc_ttt.h>
#include <init_phi.h>
#include <ode_driver.h>
#include <ire_a.h>
#include <ire_alpha.h>
/* Shapes: */
int Nx;
int shape[1];
int dim;
int level;
/* Coordinates: */
double *x;
/* Grid Functions: */
double *Trr;
double *Ttt;
double *a;
double *alpha;
double *n_phi;
double *n_pi;
double *n_pp;
double *nm1_pi;
double *nm1_pp;
double *np1_pi;
double *np1_pp;
double *lna;
double *lnalpha;
double ire_a,ire_alpha;
/* Parameters: */
double A;
double ht;
double x0;
double delx;
double epsdis;
double myzero;
double idsignum;
double totres;
double ire_pp,ire_pi;
double pp_res,pi_res;
/* Coordinate Parameters: */
double x_max;
double x_min;
double hx;
double bbox[2];
int phys_bdy[4];
/* Time Evolution Parameters: */
int steps;
int output_freq;
double lambda;
double time;

void set_cons(double *f, int N, double c) {
	int i;
	for(i=0;i<N;i++) f[i] = c;
}

void swap_levels(double **a, double **b) {
     double *t;
     t = *a;
     *a = *b;
     *b = t;
}

void copy_lev(double **a, double **b) {
	*b = *a;
}


void read_params(char *p_file,double *delx,int *Nx,int *level,double *x_max,double *x_min,double *epsdis,double *idsignum,double *lambda,double *x0,double *A,int *output_freq,int *steps,double *myzero)
{
get_param(p_file,"delx","double",1,delx);
get_param(p_file,"Nx","long",1,Nx);
get_param(p_file,"level","long",1,level);
get_param(p_file,"x_max","double",1,x_max);
get_param(p_file,"x_min","double",1,x_min);
get_param(p_file,"epsdis","double",1,epsdis);
get_param(p_file,"idsignum","double",1,idsignum);
get_param(p_file,"lambda","double",1,lambda);
get_param(p_file,"x0","double",1,x0);
get_param(p_file,"A","double",1,A);
get_param(p_file,"output_freq","long",1,output_freq);
get_param(p_file,"steps","long",1,steps);
get_param(p_file,"myzero","double",1,myzero);
}

int main(int argc, char **argv) {
char pfile[64];
int j;
strcpy(pfile,argv[1]);
/* Initialization of Coordinate: */
dim =1;
read_params(pfile,&delx,&Nx,&level,&x_max,&x_min,&epsdis,&idsignum,&lambda,&x0,&A,&output_freq,&steps,&myzero);

Nx = Nx*(int)pow(2.0,(double)level)+1;
steps = steps*(int)pow(2.0,(double)level);
/* Allocating Memory to Grid Functions: */
x = vec_alloc(Nx);
Trr = vec_alloc(1*Nx);
Ttt = vec_alloc(1*Nx);
a = vec_alloc(1*Nx);
lna = vec_alloc(1*Nx);
lnalpha = vec_alloc(1*Nx);
alpha = vec_alloc(1*Nx);
n_phi = vec_alloc(1*Nx);
n_pi = vec_alloc(1*Nx);
n_pp = vec_alloc(1*Nx);
nm1_pi = vec_alloc(1*Nx);
nm1_pp = vec_alloc(1*Nx);
np1_pi = vec_alloc(1*Nx);
np1_pp = vec_alloc(1*Nx);
hx = (x_max-x_min)/(Nx-1);
dvumsh(x,Nx,x_min,x_max);
ht = lambda*sqrt( 0.0+ hx*hx);
shape[0]=Nx;
bbox[0]=x_min;
bbox[1]=x_max;
time=0.0;

set_cons(a,Nx,1.0);
set_cons(alpha,Nx,1.0);


init_phi_(x,&Nx,&A,&x0,&delx,n_phi);

init_pi_(x,n_phi,&Nx,&hx,&myzero,&idsignum,phys_bdy,n_pi);
init_pi_(x,n_phi,&Nx,&hx,&myzero,&idsignum,phys_bdy,np1_pi);
init_pp_(x,n_phi,&Nx,&hx,&myzero,phys_bdy,n_pp);
init_pp_(x,n_phi,&Nx,&hx,&myzero,phys_bdy,np1_pp);

calc_trr_(n_pi,n_pp,&Nx,phys_bdy,Trr);
calc_ttt_(n_pi,n_pp,&Nx,phys_bdy,Ttt);

ode_driver_(lna,lnalpha,a,alpha,Ttt,Trr,&Nx,&hx,x);

gft_out_bbox("pp.sdf",time,shape,dim,bbox,n_pp);
gft_out_bbox("pi.sdf",time,shape,dim,bbox,n_pi);
gft_out_bbox("phi.sdf",time,shape,dim,bbox,n_phi);
gft_out_bbox("a.sdf",time,shape,dim,bbox,a);
gft_out_bbox("alpha.sdf",time,shape,dim,bbox,alpha);
gft_out_bbox("Trr.sdf",time,shape,dim,bbox,Trr);
gft_out_bbox("Ttt.sdf",time,shape,dim,bbox,Ttt);

ire_a_(Ttt,a,x,&Nx,&hx,&ire_a);
ire_alpha_(Trr,alpha,a,x,&Nx,&hx,&ire_alpha);
printf("step: %d t: %f iter: %d res: %1.14e ire_geom: %1.14e ire_matter: %1.14e\n",0,time,0,0.0,ire_a+ire_alpha,0.0);
ire_pp=0;
ire_pi=0;



for (int i=0; i<steps; i++) {
	j = 0;
	totres=1;
	while (totres > 1.0e-9) {
        u_pp_(a,x,alpha,n_pi,n_pp,np1_pi,np1_pp,&Nx,&ht,&hx,&epsdis,&myzero,phys_bdy,np1_pp);
        u_pi_(a,x,alpha,n_pi,n_pp,np1_pi,np1_pp,&Nx,&ht,&hx,&epsdis,phys_bdy,np1_pi);
        res_pp_(a,x,alpha,n_pi,n_pp,np1_pi,np1_pp,&Nx,&ht,&hx,&epsdis,&myzero,phys_bdy,&pp_res);
        res_pi_(a,x,alpha,n_pi,n_pp,np1_pi,np1_pp,&Nx,&ht,&hx,&epsdis,phys_bdy,&pi_res);
        totres=pp_res/l2norm(Nx,np1_pp) + pi_res/l2norm(Nx,np1_pi);
        calc_trr_(np1_pi,np1_pp,&Nx,phys_bdy,Trr);
        calc_ttt_(np1_pi,np1_pp,&Nx,phys_bdy,Ttt);
        ode_driver_(lna,lnalpha,a,alpha,Ttt,Trr,&Nx,&hx,x);
		  j++;
	}

	time = time + ht;

   ire_a_(Ttt,a,x,&Nx,&hx,&ire_a);
   ire_alpha_(a,Trr,alpha,x,&Nx,&hx,&ire_alpha);

	if ( i > 1 ) {
      ire_pi_(a,n_pp,nm1_pi,np1_pi,alpha,x,&Nx,&ht,&hx,&ire_pi);
      ire_pp_(n_pi,nm1_pp,a,alpha,np1_pp,&Nx,&ht,&hx,&ire_pp);
	}

   rdvcpy(nm1_pi,n_pi,Nx);
	rdvcpy(nm1_pp,n_pp,Nx);


	swap_levels(&np1_pp,&n_pp);
	swap_levels(&np1_pi,&n_pi);

	if ( (i + 1) % (output_freq*(int)pow(2.0,(double)level))  == 0) {
		gft_out_bbox("pp.sdf",time,shape,dim,bbox,np1_pp);
		gft_out_bbox("pi.sdf",time,shape,dim,bbox,np1_pi);
		gft_out_bbox("Trr.sdf",time,shape,dim,bbox,Trr);
		gft_out_bbox("Ttt.sdf",time,shape,dim,bbox,Ttt);
		gft_out_bbox("a.sdf",time,shape,dim,bbox,a);
		gft_out_bbox("alpha.sdf",time,shape,dim,bbox,alpha);
		printf("step: %d t: %f iter: %d res: %1.14e ire_geom: %1.14e ire_matter: %1.14e \n",i+1,time,j,totres,ire_a+ire_alpha,ire_pi+ire_pp);
	}



}
gft_close_all();
}
