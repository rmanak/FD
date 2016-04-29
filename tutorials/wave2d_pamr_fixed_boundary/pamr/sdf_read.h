#include <bbhutil.h>
#include <dv.h>
#include <iv.h>

#define max_fname_size 256
#define max_name_size 64
#define max_dim 3

struct my_SDF {
	char file[max_fname_size];
	char name[max_name_size];
	char **cnames;
	int dim;
	int size;
	int coord_size;
	int shape[max_dim];
	double bbox[max_dim*2];
	double *data;
	double *coords;
	double *x;
	double *y;
	double *z;
	double dx;
	double dy;
	double dz;
	int Nx;
	int Ny;
	int Nz;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
} ;

typedef struct my_SDF SDF;

int sdf_read(SDF *g, const char *fname, int lev, double *time);
