#include "FemTech.h"
#include "blas.h"

// plane strain or three-dimensional Linear Elastic Material
// Evaluates the PK2 stress tensor

void HourglassStress(int e, int gp) {
		double mu = properties[MAXMATPARAMS * pid[e] + 1];
		double lambda = properties[MAXMATPARAMS * pid[e] + 2];
		double E,Cshell,h,area, SRx, SRy;
		for(int i=0; i<nshell; i++)
			if(ShellID[i]==e)
				{
					h = Thickness[i];
					area = areashell[i];    		//better way to do this?
					SRx = hg_strainrate[i*nshell + 0];
					SRy = hg_strainrate[i*nshell + 1];
					break;
				}
		E = mu*(3*lambda+2*mu)/(lambda+mu);
		Cshell = 0.05*E*h*area*(Bshell[0]*Bshell[0]+Bshell[1]*Bshell[1]+Bshell[2]*Bshell[2]+Bshell[3]*Bshell[3]+Bshell[4]*Bshell[4]+Bshell[5]*Bshell[5]+Bshell[6]*Bshell[6]+Bshell[7]*Bshell[7])/8;
		for(int i=0; i<2; i++){
			HGshell_prev[HGshellptr[e]+2*gp+i]=HGshell[HGshellptr[e]+2*gp+i];
		}
		HGshell[HGshellptr[e]+2*gp+0] = HGshell_prev[HGshellptr[e]+2*gp+0] + dt*Cshell*SRx;
		HGshell[HGshellptr[e]+2*gp+1] = HGshell_prev[HGshellptr[e]+2*gp+1] + dt*Cshell*SRy;
		printf("%f %f\n", HGshell[HGshellptr[e]+2*gp+0], dt);
	return;
}
