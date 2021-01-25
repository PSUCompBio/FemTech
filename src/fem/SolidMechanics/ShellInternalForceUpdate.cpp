#include "FemTech.h"
#include "blas.h"

void ShellInternalForceUpdate(int e, int gp, double *force) {
  int nNodesL = nShapeFunctions[e];
  int bColSize = nNodesL*ndim;
  double h;
  for(int i=0; i<nshell; i++)
    {if(ShellID[i]==e)
      {
        h = Thickness[i];     //better way to do this?
        break;
      }
    }
  double *f = (double*)calloc(8, sizeof(double));
  for(int i=0; i<4; i++){
    f[i+4*0]=0.5*h*(Bshell[e*8+i+4*0]*cauchyshell[cauchyshellptr[e]+3*gp+0]+Bshell[e*8+i+4*1]*cauchyshell[cauchyshellptr[e]+3*gp+2]);
    f[i+4*1]=0.5*h*(Bshell[e*8+i+4*1]*cauchyshell[cauchyshellptr[e]+3*gp+1]+Bshell[e*8+i+4*0]*cauchyshell[cauchyshellptr[e]+3*gp+2]);
    force[i*ndim+0] = corotationalx[i*ndim]*f[i+4*0] + corotationaly[i*ndim]*f[i+4*1];
    force[i*ndim+1] = corotationalx[i*ndim + 1]*f[i+4*0] + corotationaly[i*ndim + 1]*f[i+4*1];
    force[i*ndim+2] = corotationalx[i*ndim + 2]*f[i+4*0] + corotationaly[i*ndim + 2]*f[i+4*1];
  }
	return;
  free(f);
}
