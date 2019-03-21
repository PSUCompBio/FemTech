#include "FemTech.h"
#include "blas.h"

void GetForce_3D(){
	for (int i=0;i<nelements;i++){
		//for (j = 0; j < nnel_normal; j++) {
		// g = (*elements_host)(i, j + 2);
		// x_store[i][j] = (*nodes_host)(g, 1);
		// y_store[i][j] = (*nodes_host)(g, 2);
		// z_store[i][j] = (*nodes_host)(g, 3);
	 //}
		//printf("number of GaussPoints %d\n", GaussPoints[i]);
		//printf("number of nShapeFunctions %d\n",nShapeFunctions[i]);
		//for (int ix = 0; ix < 2; ix++) {
		//	 x = points_normal(ix);
//			 wtx_normal[ix] = weights_normal(ix);

//			 for (iy = 0; iy < 2; iy++) {
//					 y = points_normal(iy);
//					 wty_normal[ix][iy] = weights_normal(iy);

//					 for (iz = 0; iz < 2; iz++) {
//							 z = points_normal(iz);
//							 wtz_normal[ix][iy][iz] = weights_normal(iz);

//							 fe_dniso_8_array(dndr_store[i][ix][iy][iz], dnds_store[i][ix][iy][iz], dndt_store[i][ix][iy][iz], x, y, z, ix, iy, iz);

//							 fe_calJacobian_array(jacobian_store[i][ix][iy][iz], nnel_normal, dndr_store[i][ix][iy][iz], dnds_store[i][ix][iy][iz], dndt_store[i][ix][iy][iz], x_store[i], y_store[i], z_store[i]);

//							 det_store[i][ix][iy][iz] = fe_detMatrix_pbr_array(jacobian_store[i][ix][iy][iz]);

//							 fe_invMatrix_pbr_array(invJacobian_store[i][ix][iy][iz], jacobian_store[i][ix][iy][iz], det_store[i][ix][iy][iz]);

//							 fe_dndxyz_8_pbr_array(dndx_store[i][ix][iy][iz], dndy_store[i][ix][iy][iz], dndz_store[i][ix][iy][iz], nnel_normal, dndr_store[i][ix][iy][iz], dnds_store[i][ix][iy][iz], dndt_store[i][ix][iy][iz], invJacobian_store[i][ix][iy][iz]);
//					 } // loop on iz
//			 } // loop on iy
	 //} //loop on ix

		//double *FT = (double*)malloc(GaussPoints[i] * ndim * sizeof(double));// transpose of F
		for(int j=0;j<GaussPoints[i];j++){
			CalculateDeformationGradient(i,j);
			DeterminateF(i,j);
			//Inverse(F,i,j);
			StressDisplacementMatrix(i,j);
			StressUpdate(i,j);
		}
    //free(FT);
	}// loop on i, nelements


	return;
}
