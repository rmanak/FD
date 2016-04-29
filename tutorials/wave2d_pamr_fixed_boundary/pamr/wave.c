//=============================================================================
// application interface functions for wave example
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <pamr.h>
#include <amrd.h>
#include <math.h>
#include <mpi.h>
#include "wave.h"

#include <sys/timeb.h>
void elapsed_time(void);

//=============================================================================
// id parameters (time symmetric gaussian)
//=============================================================================

real A, xc, yc, delx,dely,idsigx,idsigy,res_f,res_f_t;
real ev_tol;

//=============================================================================
// some convenient, "local" global variables
//=============================================================================

real *n_f,*np1_f; // in AMRH n/np1/nm1
real *n_f_t,*np1_f_t;

real hx,hy,ht;
real myzero;
real *x,*y;
int my_rank;
int shape[2],ghost_width[4],Nx,Ny,phys_bdy[4],size,g_rank,dim;
real base_bbox[4],bbox[4],dx,dy,dt;
int g_L;

int n_f_gfn, np1_f_gfn, n_f_t_gfn, np1_f_t_gfn;

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((n_f_gfn   = PAMR_get_gfn("f",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((n_f_t_gfn   = PAMR_get_gfn("f_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
	 if ((np1_f_gfn = PAMR_get_gfn("f",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((np1_f_t_gfn = PAMR_get_gfn("f_t",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
}

//=============================================================================
// call with valid iter to set up globals:
//=============================================================================
void ldptr(void)
{
   real dx0[2];
   real *x0[2],*gfs[PAMR_MAX_GFNS];
   static int first=1;

   if (first) 
   {
      first=0; 
      set_gfns();
      PAMR_get_global_bbox(base_bbox);
   }

   PAMR_get_g_dim(&dim);
   PAMR_get_g_rank(&g_rank);
   PAMR_get_g_shape(shape);
   PAMR_get_g_bbox(bbox);
   PAMR_get_g_ghost_width(ghost_width);
   PAMR_get_g_level(&g_L);
   PAMR_get_dxdt(g_L,dx0,&dt);
   ht = dt;
   dx=dx0[0];
   hx = dx;
   dy=dx0[1];
   hy = dy;

   if ((bbox[0]-base_bbox[0])<dx/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dx/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   Nx=shape[0];
   size=Nx;
   if ((bbox[2]-base_bbox[2])<dy/2) phys_bdy[2]=1; else phys_bdy[2]=0;
   if ((base_bbox[3]-bbox[3])<dy/2) phys_bdy[3]=1; else phys_bdy[3]=0;
   Ny=shape[1];
   size*=Ny;
  PAMR_get_g_x(x0);

   x=x0[0];
   y=x0[1]; 

   PAMR_get_g_gfs(gfs);

   n_f   = gfs[n_f_gfn-1];
   n_f_t   = gfs[n_f_t_gfn-1];
   np1_f = gfs[np1_f_gfn-1];
   np1_f_t = gfs[np1_f_t_gfn-1];

}

//=============================================================================
// utility routines
//=============================================================================
void const_f(real *f, real c)
{
   int i;

   for (i=0; i<Nx*Ny; i++) f[i]=c;
}

void zero(real *f)
{
   const_f(f,0);
}

//=============================================================================
// Routines required by amrd:
//=============================================================================

//=============================================================================
// Returns 0 to use default mechanism, or is expected to calculate
// the correct initial hierarchy and return 1:
//=============================================================================
int wave_id(void)
{
	if( my_rank == 0 ) elapsed_time();
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void wave_var_pre_init(char *pfile)
{
   return;
}

void wave_var_post_init(char *pfile)
{
   int i,j;
   char buf[64];

   if (my_rank==0)
   {
		system("date > date.dat");
      printf("===================================================================\n");
      printf("Reading wave parameters:\n\n");
   }


   AMRD_real_param(pfile,"A",&A,1);
   AMRD_real_param(pfile,"xc",&xc,1);
   AMRD_real_param(pfile,"yc",&yc,1);
   AMRD_real_param(pfile,"delx",&delx,1);
   AMRD_real_param(pfile,"dely",&dely,1);
   AMRD_real_param(pfile,"idsigx",&idsigx,1);
   AMRD_real_param(pfile,"idsigy",&idsigy,1);
   AMRD_real_param(pfile,"myzero",&myzero,1);
   AMRD_real_param(pfile,"ev_tol",&ev_tol,1);

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all t=n variables to their 'zero' values:
//=============================================================================
void wave_AMRH_var_clear(void)
{
   ldptr();

   zero(n_f);  zero(n_f_t);

   return;
}

//=============================================================================
// Initial data for free fields: (at tn=2)
//
// currently, we only allow for time-symmetric initial data
//=============================================================================
void wave_free_data(void)
{
   ldptr();

	// Initializer:
   init_f_(y,x,&Nx,&Ny,&A,&xc,&yc,&delx,&dely,n_f);
   init_f_(y,x,&Nx,&Ny,&A,&xc,&yc,&delx,&dely,np1_f);
   init_f_t_(y,x,&Nx,&Ny,&A,&xc,&yc,&delx,&dely,&idsigx,&idsigy,n_f_t);
   init_f_t_(y,x,&Nx,&Ny,&A,&xc,&yc,&delx,&dely,&idsigx,&idsigy,np1_f_t);

   return;
}  

//=============================================================================
// Initial constraint data --- called after each MG iteration.
//=============================================================================
void wave_t0_cnst_data(void)
{
   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//
// NOTE: at this point, the time sequence is: n,nm1,np1
//=============================================================================
void wave_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration.
// We're using an explicit scheme to solve for phi, hence return 0
//=============================================================================
real wave_evo_residual(void)
{
   // Residual:
	ldptr();

   res_f_(n_f,x,y,n_f_t,np1_f,np1_f_t,&Nx,&Ny,&ht,&myzero,phys_bdy,&res_f);

   res_f_t_(n_f,x,y,n_f_t,np1_f,np1_f_t,&Nx,&Ny,&ht,&hx,&hy,&myzero,phys_bdy,&res_f_t);

   return res_f+res_f_t;
}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
real wave_MG_residual(void)
{
   return 0;
}

//=============================================================================
// Performs 1 iteration of the evolution equations 
//=============================================================================
void wave_evolve(int iter)
{
   ldptr();
   double tot_res;

// Evolution 1 step:
    
    tot_res = 1.0;
    while ( tot_res > ev_tol ) { 
    u_f_(n_f,x,y,n_f_t,np1_f,np1_f_t,&Nx,&Ny,&ht,&myzero,phys_bdy,np1_f);
    u_f_t_(n_f,x,y,n_f_t,np1_f,np1_f_t,&Nx,&Ny,&ht,&hx,&hy,&myzero,phys_bdy,np1_f_t);
    res_f_(n_f,x,y,n_f_t,np1_f,np1_f_t,&Nx,&Ny,&ht,&myzero,phys_bdy,&res_f);
    res_f_t_(n_f,x,y,n_f_t,np1_f,np1_f_t,&Nx,&Ny,&ht,&hx,&hy,&myzero,phys_bdy,&res_f_t);
    tot_res = res_f_t+res_f;
    }

   return;

}

//=============================================================================
// sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS!!!)
//=============================================================================
void wave_fill_ex_mask(real *mask, int dim, int *shape, real *bbox, real excised)
{
}

//=============================================================================
void wave_fill_bh_bboxes(real *bbox, int *num, int max_num)
{
}

//=============================================================================
void wave_post_tstep(int L)
{
   return;
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
real wave_MG_relax(void)
{
   return 0;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void wave_L_op(void)
{
   return;
}

//=============================================================================
// Called after calculating the TRE for all variables
//=============================================================================
void wave_scale_tre(void)
{
   return;
}

//=============================================================================
// post-regrid initialization of constant functions
//=============================================================================
void wave_post_regrid(void)
{
}

//=============================================================================
int main(int argc, char **argv)
{
   amrd(argc,argv,&wave_id,&wave_var_pre_init,
        &wave_var_post_init, &wave_AMRH_var_clear,
        &wave_free_data, &wave_t0_cnst_data,
        &wave_evo_residual, &wave_MG_residual,
        &wave_evolve, &wave_MG_relax, &wave_L_op, 
        &wave_pre_io_calc, &wave_scale_tre, 
        &wave_post_regrid, &wave_post_tstep,
        &wave_fill_ex_mask, &wave_fill_bh_bboxes);
   if (my_rank==0) elapsed_time();
}

//=============================================================================
// Maintains/reports elapsed wall-clock time.
//=============================================================================
void elapsed_time(void) {
   static int    first = 1;
   struct        timeb t;
   static double msinit;
   double        mscurr, mselapsed;

   ftime(&t);
   mscurr = 1000.0 * t.time + t.millitm;
   if( first ) {
      msinit = mscurr;
      first = 0;
   }
	mselapsed = mscurr - msinit;
   printf("elapsed_time: Seconds since initial call: %12.3f\n",
         mselapsed / 1000.0);
}
//ire_f_(np1_f,n_f,nm1_f,&Nx,&Ny,&ht,&hx,&hy,res);
