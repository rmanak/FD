#ifndef _NUM
#define _NUM 1
void init_f_(double *y,double *x,int *Nx,int *Ny,double *A,double *xc,double *yc,double *delx,double *dely,double *res);

void init_f_t_(double *y,double *x,int *Nx,int *Ny,double *A,double *xc,double *yc,double *delx,double *dely,double *idsigx,double *idsigy,double *res);

void ire_f_(double *np1_f,double *n_f,double *nm1_f,int *Nx,int *Ny,double *ht,double *hx,double *hy,double *res);

void res_f_(double *n_f,double *x,double *y,double *n_f_t,double *np1_f,double *np1_f_t,int *Nx,int *Ny,double *ht,double *myzero,int *phys_bdy,double *res);

void res_f_t_(double *n_f,double *x,double *y,double *n_f_t,double *np1_f,double *np1_f_t,int *Nx,int *Ny,double *ht,double *hx,double *hy,double *myzero,int *phys_bdy,double *res);

void u_f_(double *n_f,double *x,double *y,double *n_f_t,double *np1_f,double *np1_f_t,int *Nx,int *Ny,double *ht,double *myzero,int *phys_bdy,double *res);

void u_f_t_(double *n_f,double *x,double *y,double *n_f_t,double *np1_f,double *np1_f_t,int *Nx,int *Ny,double *ht,double *hx,double *hy,double *myzero,int *phys_bdy,double *res);

#endif
