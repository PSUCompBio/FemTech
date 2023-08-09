#include "FemTech.h"
#include "blas.h"
#include "lapack.h"

int fibermatrixID(){
/*	double elementCoordinates[24], felementCoordinates[6];
        int ele;	
	for (int i=0; i<nembedel; i++){
		felid = embedel[i];
		for(int j=0; j<3; j++){
			felementCoordinates[j] = coordinates[ndim*connectivity[eptr[felid]]+j];
			felementCoordinates[j+3] = coordinates[ndim*connectivity[eptr[felid+1]]+j];
			}
		for(int j=0; j<nelements-nembedel; j++){
			int elid = hostel[j];
			for (int k = 0; k < 3; ++k) {
      				int index = ndim*connectivity[eptr[elid]]+k;
      				elementCoordinates[0*ndim+k] = coordinates[index];
  			}		
	const int index[24] = {0, 3, 2, 1, 4, 5, 6, 7, 0, 4, 7, 3, 1, 2, 6, 5, 0, 1, 5, 4, 3, 7, 6, 2};
  for(int m=0; m<6; m++){
  double p1[3], p2[3], p3[3], p4[3];
  p1[0] = (elementCoordinates[3*index[m*4+0]]);
  p2[0] = (elementCoordinates[3*index[m*4+1]]);
  p3[0] = (elementCoordinates[3*index[m*4+2]]);
  p4[0] = (elementCoordinates[3*index[m*4+3]]);
  p1[1] = (elementCoordinates[3*index[m*4+0]+1]);
  p2[1] = (elementCoordinates[3*index[m*4+1]+1]);
  p3[1] = (elementCoordinates[3*index[m*4+2]+1]);
  p4[1] = (elementCoordinates[3*index[m*4+3]+1]);
  p1[2] = (elementCoordinates[3*index[m*4+0]+2]);
  p2[2] = (elementCoordinates[3*index[m*4+1]+2]);
  p3[2] = (elementCoordinates[3*index[m*4+2]+2]);
  p4[2] = (elementCoordinates[3*index[m*4+3]+2]);
  double dist1[3], dist2[3], dist3[3], dist4[3], normal1[3], normal2[3];
  for (int i = 0; i < 3; ++i) {
        dist1[i] = (p2[i]-p1[i]);
        dist2[i] = (p3[i]-p2[i]);
        dist3[i] = (p4[i]-p3[i]);
	dist4[i] = (p4[i]-p1[i]);
        }
  crossProduct(dist1, dist2, normal1);
  crossProduct(dist3, dist4, normal2);
  }
  }
}
  return ele;*/
}
