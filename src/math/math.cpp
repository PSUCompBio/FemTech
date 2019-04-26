#include "FemTech.h"
#include "blas.h"

int math()
{
   // printf() displays the string inside quotation
   printf("Hello, World and Digital Brain!!\n");
   return 0;
}

void inverse3x3Matrix(double* mat, double* invMat, double* det) {
  // Mat and invMat are 1d arrays with colum major format for storing matrix
  // Compute matrix determinant
  double detLocal = mat[0] * (mat[4] * mat[8] - mat[5] * mat[7]) -
              mat[3] * (mat[1] * mat[8] - mat[7] * mat[2]) +
              mat[6] * (mat[1] * mat[5] - mat[4] * mat[2]);

  double invdet = 1 / detLocal;

  invMat[0] = (mat[4] * mat[8] - mat[5] * mat[7]) * invdet;
  invMat[3] = (mat[6] * mat[5] - mat[3] * mat[8]) * invdet;
  invMat[6] = (mat[3] * mat[7] - mat[6] * mat[4]) * invdet;
  invMat[1] = (mat[7] * mat[2] - mat[1] * mat[8]) * invdet;
  invMat[4] = (mat[0] * mat[8] - mat[6] * mat[2]) * invdet;
  invMat[7] = (mat[1] * mat[6] - mat[0] * mat[7]) * invdet;
  invMat[2] = (mat[1] * mat[5] - mat[2] * mat[4]) * invdet;
  invMat[5] = (mat[2] * mat[3] - mat[0] * mat[5]) * invdet;
  invMat[8] = (mat[0] * mat[4] - mat[1] * mat[3]) * invdet;

  (*det) = detLocal;
}

double normOfCrossProduct(double *a, double *b) {
  double z = a[0]*b[1]-a[1]*b[0];
  if (ndim == 3) {
    double x = a[1]*b[2]-a[2]*b[1];
    double y = -a[0]*b[2]+a[2]*b[0];
    return sqrt(x*x + y*y + z*z);
  } else {
    if (ndim == 2) {
      return z;
    }
  }
  return 0.0;
}

double tripleProduct(double *s, double *a, double *b) {
  return s[2]*(a[0]*b[1]-a[1]*b[0])+
         s[0]*(a[1]*b[2]-a[2]*b[1])-
         s[1]*(a[0]*b[2]-a[2]*b[0]);
}

void crossProduct(double* a, double* b, double* result) {
  result[0] = a[1]*b[2]-a[2]*b[1];
  result[1] = -a[0]*b[2]+a[2]*b[0];
  result[2] = a[0]*b[1]-a[1]*b[0];
}
double norm3D(double *a) {
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}
double dotProduct3D(double *a, double *b) {
  return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}
/*
void MultiplyMatrices(double* a[], double* b[], int sizeM, double* result ){
	//Multiplication Logic
	int sum;
	double c[sizeM][sizeM];
	for (int i = 0; i < sizeM; i++) {
		 for (int j = 0; j < sizeM; j++) {
				sum = 0;
				for (int k = 0; k < sizeM; k++) {
					 sum = sum + a[i][k] * b[k][j];
				}
				c[i][j] = sum;
		 }
	}

	(*result) = c;
}
*/
