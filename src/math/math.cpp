#include "FemTech.h"
#include "blas.h"
#include "lapack.h"

#include <algorithm>
#include <vector>
#include <cstdlib>

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
// n and x are assumed to be three dimensional
// Source: http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf
// After rotation the results are stored in the input matrix
// Theta is assumed to be in radians
void rotate3d(double *n, double theta, double *xin) {
  double x, y, z;
  const double cth = cos(theta);
  const double sth = sin(theta);
  const double mcth_m1 = 1.0-cth;
  x = (cth+n[0]*n[0]*mcth_m1)*xin[0] + (n[0]*n[1]*mcth_m1-n[2]*sth)*xin[1] + \
      (n[0]*n[2]*mcth_m1+n[1]*sth)*xin[2];
  y = (cth+n[1]*n[1]*mcth_m1)*xin[1] + (n[1]*n[2]*mcth_m1-n[0]*sth)*xin[2] + \
      (n[0]*n[1]*mcth_m1+n[2]*sth)*xin[0];
  z = (cth+n[2]*n[2]*mcth_m1)*xin[2] + (n[0]*n[2]*mcth_m1-n[1]*sth)*xin[0] + \
      (n[2]*n[1]*mcth_m1+n[0]*sth)*xin[1];
  xin[0] = x; xin[1] = y; xin[2] = z;
}
// Create the rotation matrix in row major format for faster multiplication with
// position vector
void get3dRotationMatrix(double *n, double theta, double mat[3][3]) {
  const double cth = cos(theta);
  const double sth = sin(theta);
  const double mcth_m1 = 1.0-cth;

  mat[0][0] = cth+n[0]*n[0]*mcth_m1;
  mat[0][1] = n[0]*n[1]*mcth_m1-n[2]*sth;
  mat[0][2] = n[0]*n[2]*mcth_m1+n[1]*sth;

  mat[1][0] = n[0]*n[1]*mcth_m1+n[2]*sth;
  mat[1][1] = cth+n[1]*n[1]*mcth_m1;
  mat[1][2] = n[1]*n[2]*mcth_m1-n[0]*sth;

  mat[2][0] = n[0]*n[2]*mcth_m1-n[1]*sth;
  mat[2][1] = n[2]*n[1]*mcth_m1+n[0]*sth;
  mat[2][2] = cth+n[2]*n[2]*mcth_m1;
}

// Linear interpolation of tabulated data
double interpolateLinear(int n, double *x, double *y, double value) {
  if (value < x[0]) {
    FILE_LOG_SINGLE(ERROR, "Linear Interpolation value outside range");
    return 0.0;
  }
  if (value > x[n-1]) {
    return y[n-1];
  }
  if (value == x[0]) {
    return y[0];
  }
  int index = 0;
  for (int i = 1; i < n; ++i) {
    if (value <= x[i]) {
      index = i-1;
      break;
    }
  }
  const double yValue = y[index] + (y[index+1]-y[index])*(value-x[index])/(x[index+1]-x[index]);
  return yValue;
}

// Functions for using quaternions
void quaternionExp(double *q1, double *q2) {
  double vMag = q1[1]*q1[1]+q1[2]*q1[2]+q1[3]*q1[3];
  if (vMag == 0) {
    q2[0] = 1.0; q2[1] = 0.0; q2[2] = 0.0; q2[3] = 0.0;
  } else {
    vMag = sqrt(vMag);
    const double d1 = exp(q1[0]);
    const double d2 = d1*sin(vMag)/vMag;
    q2[0] = d1*cos(vMag); q2[1] = d2*q1[1]; q2[2] = d2*q1[2]; q2[3] = d2*q1[3];
  }
}

void quaternionMultiply(double *q1, double *q2, double *qr) {
  qr[0] = q1[0]*q2[0]-q1[1]*q2[1]-q1[2]*q2[2]-q1[3]*q2[3];
  qr[1] = q1[0]*q2[1]+q1[1]*q2[0]+q1[2]*q2[3]-q1[3]*q2[2];
  qr[2] = q1[0]*q2[2]-q1[1]*q2[3]+q1[2]*q2[0]+q1[3]*q2[1];
  qr[3] = q1[0]*q2[3]+q1[1]*q2[2]-q1[2]*q2[1]+q1[3]*q2[0];
}

void quaternionInverse(double *q, double *qinv) {
  double norm = q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3];
  qinv[0] = q[0]/norm; qinv[1] = -q[1]/norm; 
  qinv[2] = -q[2]/norm; qinv[3] = -q[3]/norm;
}

void quaternionRotate(double *v, double *R, double* vp) {
  double Rv[4], Rinv[4];
  quaternionMultiply(R, v, Rv);
  quaternionInverse(R, Rinv);
  quaternionMultiply(Rv, Rinv, vp);
}

void quaternionRotate(double *v, double *R, double *Rinv, double* vp) {
  double Rv[4];
  quaternionMultiply(R, v, Rv);
  quaternionMultiply(Rv, Rinv, vp);
}

double compute95thPercentileValueBruteForce(double* dataArray, int localSize) {
  // Compute total size
  int totalSize = 0;
  double *fullData;
  int *recvCount, *recvDisplacement;
  // Calculate the total size
  MPI_Reduce(&localSize, &totalSize, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (world_rank == 0) {
    fullData = (double*)malloc(totalSize*sizeof(double));
    recvCount = (int*)malloc(world_size*sizeof(int));
    recvDisplacement = (int*)malloc(world_size*sizeof(int));
  }
  // Get recv count
  MPI_Gather(&localSize, 1, MPI_INT, recvCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (world_rank == 0) {
    recvDisplacement[0] = 0;
    for (int i = 1; i < world_size; ++i) {
      recvDisplacement[i] = recvCount[i-1] + recvDisplacement[i-1];
    }
  }
  // Get full array
  MPI_Gatherv(dataArray, localSize, MPI_DOUBLE, fullData, recvCount, \
      recvDisplacement, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  double nthElement = 0.0;
  if (world_rank == 0) {
    int index95 = (int)(totalSize*0.95)-1; //Runtime error if index95 < 0
    std::vector<double> fullDataVec(fullData, fullData + totalSize);

    std::nth_element(fullDataVec.begin(), fullDataVec.begin()+index95, fullDataVec.end());
    nthElement = fullDataVec[index95];

    free(fullData);
    free(recvCount);
    free(recvDisplacement);
  }
  MPI_Bcast(&nthElement, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return nthElement;
}

double computeMedianLocal(double *data, int size) {
  double median;
  // convert array to vector
  std::vector<double> dataVec(data, data + size);
  int mid = size/2;
  if (size%2) {
    // array size even
    std::nth_element(dataVec.begin(), dataVec.begin()+mid, dataVec.end());
    std::nth_element(dataVec.begin(), dataVec.begin()+(mid-1), dataVec.end());
    median = 0.5*(dataVec[mid]+dataVec[mid-1]);
  } else {
    // array size odd
    std::nth_element(dataVec.begin(), dataVec.begin()+mid, dataVec.end());
    median = dataVec[mid];
  }
  return median;
}

double computeWeightedMedian(double *data, int *weight, int size) {
  std::vector< std::pair<double, int> > vecPair;
  double weightSum = 0.0;
  for (int i = 0; i < size; ++i) {
    vecPair.push_back(std::make_pair(data[i], weight[i]));
    weightSum += (double)weight[i];
  }
  std::sort(vecPair.begin(), vecPair.end());

  int j = 0;
  double partialSum = weightSum - (double)weight[0];
  while(partialSum > weightSum*0.5) {
    j = j + 1;
    partialSum = partialSum - (double)weight[j];
  }
  return data[j];
}

double compute95thPercentileValue(double* dataArray, int localSize) {
  /* Implementation is based on
  A Note on Parallel Selection on Coarse-Grained Multicomputers
  E. L. G. Saukasi and S. W. Song Algorithmica (1999) 24:  371–380 */
  // Copy for use inside loop
  double *dataArrayLoop = (double*)malloc(localSize*sizeof(double));
  double *dataArrayLoopCopy = (double*)malloc(localSize*sizeof(double));
  memcpy(dataArrayLoop, dataArray, localSize*sizeof(double));
  int localSizeLoop = localSize;
  // Compute total size
  int totalSize = 0;
  MPI_Allreduce(&localSize, &totalSize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  int k = (int)(0.95*totalSize)-1;

  double *median = NULL;
  int *sizes = NULL;
  if (world_rank == 0) {
    median = (double*)malloc(world_size*sizeof(double));
    sizes = (int*)malloc(world_size*sizeof(int));
  }

  int N = totalSize;
  const int c = 20; // TODO(Anil) tune c for NP = 16
  double value;
  while (N < totalSize/(c*world_size)) {
    // compute median of local data
    double localMedian = computeMedianLocal(dataArrayLoop, localSizeLoop);
    // All processors send its median to master
    MPI_Gather(&localMedian, 1, MPI_DOUBLE, median, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // All processors send its number of elements to master as weights
    MPI_Gather(&localSizeLoop, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    double weightedMedian;
    if (world_rank == 0) {
      // Master computes the weigted median
      weightedMedian = computeWeightedMedian(median, sizes, world_size);
    }
    // Master shares the computed weighted median to all processors
    MPI_Bcast(&weightedMedian, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int li = 0, ei = 0, gi = 0;
    int L = 0, E = 0, G = 0;
    for (int i = 0; i < localSizeLoop; ++i) {
      double element = dataArrayLoop[i];
      if (element < weightedMedian) {
        li = li + 1;
      } else {
        if (element == weightedMedian) {
          ei = ei + 1;
        } else {
          gi = gi + 1;
        }
      }
    }
    // Share li, ei and gi with master
    MPI_Allreduce(&li, &L, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ei, &E, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&gi, &G, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (L < k && k <= (L+E)) {
      value = weightedMedian;
      break;
    } else {
      if (k < L) {
        N = L;
        memcpy(dataArrayLoopCopy, dataArrayLoop, sizeof(double)*localSizeLoop);
        int count = 0;
        for (int i = 0; i < localSizeLoop; ++i) {
          if (dataArrayLoopCopy[i] < weightedMedian) {
            dataArrayLoop[count] = dataArrayLoopCopy[i];
            count = count + 1;
          }
        }
        localSizeLoop = count;
      } else {
        N = G;
        k = k - (L+E);
        memcpy(dataArrayLoopCopy, dataArrayLoop, sizeof(double)*localSizeLoop);
        int count = 0;
        for (int i = 0; i < localSizeLoop; ++i) {
          if (dataArrayLoopCopy[i] > weightedMedian) {
            dataArrayLoop[count] = dataArrayLoopCopy[i];
            count = count + 1;
          }
        }
        localSizeLoop = count;
      }
    }
  }
  // Send remaining N elements to master
  value = compute95thPercentileValueBruteForce(dataArrayLoop, localSizeLoop);

  if (world_rank == 0) {
    free(median);
    free(sizes);
  }
  free(dataArrayLoop);
  free(dataArrayLoopCopy);
  return value;
}

int compareDouble(const void *a, const void *b) { 
  double diff = (*(double *)a - *(double *)b); 
  if (diff < 0) {
    return -1;
  } else if (diff > 0) {
    return 1;
    } else {
      return 0;
    }
}

void test95Percentile(void) {
  std::srand(126);
  int baseSize = 1011;
  int increment = 20;
  int count = baseSize*world_size+(increment*(world_size-1)*world_size)/2;
  double *randArray;
  int localSize = baseSize + increment*world_rank;
  double *localArray = (double*)malloc(localSize*sizeof(double));
  int *sendSize, *sendLocation;
  if (world_rank == 0) {
    randArray = (double*)malloc(count*sizeof(double));
    for (int i = 0; i < count; ++i) {
      randArray[i] = 2.0*(double)std::rand()/(double)RAND_MAX;
    }
    sendSize = (int*)malloc(sizeof(int)*world_size);
    sendLocation = (int*)malloc(sizeof(int)*world_size);
    for (int i = 0; i < world_size; ++i) {
      sendSize[i] = baseSize + increment*i;
    }
    sendLocation[0] = 0;
    for (int i = 1; i < world_size; ++i) {
      sendLocation[i] = sendLocation[i-1]+sendSize[i-1];
    }
  }
  MPI_Scatterv(randArray, sendSize, sendLocation, MPI_DOUBLE, localArray, localSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  double computed = compute95thPercentileValueBruteForce(localArray, localSize);
  double computed2 = compute95thPercentileValue(localArray, localSize);
  double actual = 0;
  if (world_rank == 0) {
    std::qsort(randArray, count, sizeof(double), compareDouble);
    int index95 = (int)(count*0.95-1.0);
    actual = randArray[index95];
  }
  FILE_LOG(INFO, "95 Percentile Test Direct   : computed = %16.9e, actual = %16.9e", computed, actual);
  FILE_LOG(INFO, "95 Percentile Test Parallel : computed = %16.9e, actual = %16.9e", computed2, actual);

  free(localArray);
  if (world_rank == 0) {
    free(randArray);
    free(sendSize);
    free(sendLocation);
  }
}

// Computes mat = mat + \alpha*(x*x^T)
void dyadic(const double* const vector, const double preFactor, double * const mat) {
  dger_(&ndim, &ndim, &preFactor, vector, &oneI, vector, &oneI, mat, &ndim);
}
