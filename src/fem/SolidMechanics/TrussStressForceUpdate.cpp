#include "FemTech.h"
#include "blas.h"

void TrussStressForceUpdate(int e, int gp, double *force){
      // int e :: element ID
      // int gp :: Gauss point ID
      // double *force :: the truss force (internal efforts)

      //
      // initial coordinates
      //
      double x[2] = {0};
      double y[2] = {0};
      double z[2] = {0};
      //
      // displacements
      //
      double dx[2] = {0};
      double dy[2] = {0};
      double dz[2] = {0};

      for (int k=0; k < eptr[e+1]; k++){
            x[k] = coordinates[ndim*connectivity[eptr[e]+k]+0];
            y[k] = coordinates[ndim*connectivity[eptr[e]+k]+1];
            z[k] = coordinates[ndim*connectivity[eptr[e]+k]+2];

            dx[k] = x[k] + displacements[ndim*connectivity[eptr[e]+k]+0];
            dy[k] = y[k] + displacements[ndim*connectivity[eptr[e]+k]+1];
            dz[k] = z[k] + displacements[ndim*connectivity[eptr[e]+k]+2];
      }

      //
      // Turn the 6-dof to 2-dof truss
      //
      double cx = 0.0;
      double cy = 0.0;
      double cz = 0.0;

      get_cos(e, &cx, &cy, &cz);

      double dl[2] = {0};

      for (int k=0; k < eptr[e+1]; k++){
            dl[k] = dx[k]*cx + dy[k]*cy + dz[k]*cz;
      }

      //
      // Green-Lagrange strain tensor
      //
      int wIndex = gpPtr[e] + gp;
      const int indexDshp = dsptr[e] + gp * GaussPoints[e] * ndim;
      const double preFactor = gaussWeights[wIndex]*detJacobian[wIndex];

      double Du[2] = {0};
      double Eu[2] = {0}; /*Since its always diagonal store just two of them*/

      for (int k=0; k < eptr[e+1]; k++){
            Du[k] = Du[k] + dshp[indexDshp + ndim * k + 0] * preFactor * dl[k]; /*We only use the x-derivatives*/
            Eu[k] = Du[k] + 0.5 * pow(Du[k],2);
      }

      //
      // Define here the local forces (2-dofs)
      //
      double fi[2] = {0};
      double young = 1.0;
      double area = 1.0;
      for (int k=0; k < eptr[e+1]; k++){
            fi[k] = fi[k] + young * area * Eu[k];
      }

      //
      // Reconstruct here the local forces (6-dofs) Need to review this
      //
      force[0] = cx * fi[0];
      force[1] = cy * fi[0];
      force[2] = cz * fi[0];

      force[3] = cx * fi[1];
      force[4] = cy * fi[1];
      force[5] = cz * fi[1];

	if(debug && 1==0){
		printf("--------F for gauss point %d --------\n",gp);
		for(int i=0;i<ndim;i++){
			for(int j=0;j<ndim;j++){
					int index = fptr[e] + ndim*ndim*gp + ndim*j+i;
					printf("%3.3e   ",F[index]);
			} //loop on j
			printf("\n");
		} //loop on i
	} // if debug

	return ;
}
