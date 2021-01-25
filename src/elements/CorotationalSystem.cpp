#include "FemTech.h"
#include <math.h>
#include <string>

double* r31;
double* r42;
double* r21;
double* s33;
double* s11;
double length;
double r21e3;
double* e1;
double* e2;
double* e3;
const int csize = ndim * sizeof(double);

void CorotationalSystem() {
	r31=(double*)calloc(ndim,sizeof(double));
	r42=(double*)calloc(ndim,sizeof(double));
	r21=(double*)calloc(ndim,sizeof(double));
	s33=(double*)calloc(ndim,sizeof(double));
	s11=(double*)calloc(ndim,sizeof(double));
	e1=(double*)calloc(ndim,sizeof(double));
	e2=(double*)calloc(ndim,sizeof(double));
	e3=(double*)calloc(ndim,sizeof(double));
  if (!corotationalx || !corotationaly || !corotationalz) {
    FILE_LOG_SINGLE(ERROR, "Error in allocating space for corotational system");
    TerminateFemTech(12);
  }
 	for(int i=0; i<nshell; i++)
 		{
			r31[0] = coordinates[connectivity[eptr[ShellID[i]]+2]*ndim] - coordinates[connectivity[eptr[ShellID[i]]]*ndim];
			r31[1] = coordinates[connectivity[eptr[ShellID[i]]+2]*ndim+1] - coordinates[connectivity[eptr[ShellID[i]]]*ndim+1];
			r31[2] = coordinates[connectivity[eptr[ShellID[i]]+2]*ndim+2] - coordinates[connectivity[eptr[ShellID[i]]]*ndim+2];
			r42[0] = coordinates[connectivity[eptr[ShellID[i]]+3]*ndim] - coordinates[connectivity[eptr[ShellID[i]]+1]*ndim];
			r42[1] = coordinates[connectivity[eptr[ShellID[i]]+3]*ndim+1] - coordinates[connectivity[eptr[ShellID[i]]+1]*ndim+1];
			r42[2] = coordinates[connectivity[eptr[ShellID[i]]+3]*ndim+2] - coordinates[connectivity[eptr[ShellID[i]]+1]*ndim+2];
			r21[0] = coordinates[connectivity[eptr[ShellID[i]]+1]*ndim] - coordinates[connectivity[eptr[ShellID[i]]]*ndim];
			r21[1] = coordinates[connectivity[eptr[ShellID[i]]+1]*ndim+1] - coordinates[connectivity[eptr[ShellID[i]]]*ndim+1];
			r21[2] = coordinates[connectivity[eptr[ShellID[i]]+1]*ndim+2] - coordinates[connectivity[eptr[ShellID[i]]]*ndim+2];
			crossProduct(r31, r42, s33);
			length = norm3D(s33);
			e3[0] = s33[0]/length;
			e3[1] = s33[1]/length;
			e3[2] = s33[2]/length;
			r21e3 = dotProduct3D(r21, corotationalz);
			s11[0] = r21[0] - r21e3*corotationalz[i*ndim];
			s11[1] = r21[1] - r21e3*corotationalz[i*ndim+1];
			s11[2] = r21[2] - r21e3*corotationalz[i*ndim+2];
			length = norm3D(s11);
			e1[0] = s11[0]/length;
			e1[1] = s11[1]/length;
			e1[2] = s11[2]/length;
			crossProduct(e3, e1, e2);
			corotationalx[i*ndim] = e1[0];
			corotationalx[i*ndim+1] = e1[1];
			corotationalx[i*ndim+2] = e1[2];
			corotationaly[i*ndim] = e2[0];
			corotationaly[i*ndim+1] = e2[1];
			corotationaly[i*ndim+2] = e2[2];
			corotationalz[i*ndim] = e3[0];
			corotationalz[i*ndim+1] = e3[1];
			corotationalz[i*ndim+2] = e3[2];
		}
		free(r31);
		free(r42);
		free(r21);
		free(s33);
		free(s11);
		free(e1);
		free(e2);
		free(e3);
	return;
}
//coordinates[connectivity[eptr[ShellID[i]]+2]*ndim+1] :shell id is the element id for shell, eptr points to the first node in that element,connectivity lists that node id, coordinates gives the actual location
