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
