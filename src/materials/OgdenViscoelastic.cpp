#include "FemTech.h"
#include "blas.h"
#include "lapack.h"

#include <stdlib.h>

/* Implements the Ogden Model
 * Expected material properties :
 * \rho      = properties(0)
 * K         = properties(1)
 * N_{Ogden} = properties(2), maximum value of 3
 * \alpha_i  = properties(3+2*i) i = 0 to N_{Ogden}-1
 * \mu_i     = properties(4+2*i) i = 0 to N_{Ogden}-1
 * N_{Prony} = properties(3+N_{Ogden}*2), maximum value of 6
 * g_i       = properties(4+2*N_{Ogden}+2*j) j = 0 to N_{Prony}-1
 * \tau_i    = properties(5+2*N_{Ogden}+2*j) j = 0 to N_{Prony}-1
 * MAXMATPARAMS = 4+2*(N_{Ogden}+N_{Prony}) = 22
 * */

void OgdenViscoelastic(int e, int gp) {
  if (ndim == 2) {
    FILE_LOG_SINGLE(ERROR, "Plane Strain implementation yet to be done");
    TerminateFemTech(3);
  }
  // Assumes ndim == 3

  /* Read the material properties */
  // Get material property ID to read material properties
  const int pideIndex = pid[e]*MAXMATPARAMS;
  // Get location of array with material properties of the element
  const double* localProperties = &(properties[pideIndex]);
  // Read material properties
  const double K = localProperties[1];
  const int nTerm = static_cast<int>(localProperties[2]);
  const int nPronyLocal = int(localProperties[3+nTerm*2]); // No of terms in prony series
  // TODO(Anil) convert local allocations to one time global
  double *gi, *ti;
  gi = (double*)malloc(nPronyLocal*sizeof(double));
  ti = (double*)malloc(nPronyLocal*sizeof(double));
  for (unsigned int i = 0; i < nPronyLocal; ++i) {
    gi[i] = localProperties[4+2*(i+nTerm)];
    ti[i] = localProperties[5+2*(i+nTerm)];
  }

  // C = F^T F
  double *H = mat1;
  ComputeH(e, gp, H);
  // Compute and store F = H + I
  const unsigned int index = fptr[e] + ndim * ndim * gp;
  double * const F_element_gp = &(F[index]);
  for (unsigned int i = 0; i < ndim2; ++i) {
    F_element_gp[i] = H[i];
  }
  F_element_gp[0] = F_element_gp[0] + 1.0;
  F_element_gp[4] = F_element_gp[4] + 1.0;
  F_element_gp[8] = F_element_gp[8] + 1.0;
  
  double *Cmat = mat2;
  // compute C = F^T F
  dgemm_(chy, chn, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
      F_element_gp, &ndim, &zero, Cmat, &ndim);

  // Compute eigen values and eigenvectors of C
  double cEigenValue[ndim];
  double dWork[workSize];
  dsyev_(jobzV, uploU, &ndim, Cmat, &ndim, cEigenValue, dWork, &workSize, &info);
  // Cmat now contains the eigen vectors of C

  // Correct Eigen Values for numerical inaccuracy
  double delat01 = fabs(cEigenValue[0] - cEigenValue[1]);
  double delat12 = fabs(cEigenValue[1] - cEigenValue[2]);
  double delat02 = fabs(cEigenValue[0] - cEigenValue[2]);

  const double eigenTol = 1e-12;
  if (delat01 < eigenTol) cEigenValue[1] = cEigenValue[0];
  if (delat02 < eigenTol) cEigenValue[2] = cEigenValue[0];
  if (delat12 < eigenTol) cEigenValue[2] = cEigenValue[1];

  const double J = det3x3Matrix(F_element_gp);
  const double Jm13 = pow(J, -1.0/3.0);
  const double invJ = 1.0/J;

  // Compute parameters derived from the eigen value required for stress
  // computations
  double eigenPower[ndim*nTerm];
  double eigenPowerSum[nTerm];
  for (int j = 0; j < nTerm; ++j) {
    eigenPowerSum[j] = 0.0;
  }
  for (int i = 0; i < ndim; ++i) {
    // Compute sqrt of lambda and multiply by J^{-1/3}
    const double lambdaBar = Jm13*sqrt(cEigenValue[i]);
    for (int j = 0; j < nTerm; ++j) {
      const double alpha = localProperties[3+j*2];
      const double lambdaBarAlpha = pow(lambdaBar, alpha);
      eigenPowerSum[j] += lambdaBarAlpha;
      eigenPower[i*nTerm+j] = lambdaBarAlpha;
    }
  }
  for (int j = 0; j < nTerm; ++j) {
    eigenPowerSum[j] /= 3.0;
  }

  // Compute the distortional/deviatoric part of PK2 stress
  double *STilde = mat3;
  for(int i = 0; i< ndim2; ++i) {
    STilde[i] = 0.0;
  }
  for (int i = 0; i < ndim; ++i) {
    double preFactor = 0.0;
    for (int j = 0; j < nTerm; ++j) {
      double mu = localProperties[4+j*2];
      preFactor += mu*(eigenPower[i*nTerm+j]-eigenPowerSum[j]);
    }
    preFactor *= invJ;
    dyadic(&Cmat[3*i], preFactor, STilde);
  }
  FILE_LOG_SINGLE(WARNING, "Trace of STilde = %15.9f", STilde[0]+STilde[4]+STilde[8]);

  // Transform Stilde to sigmaTilde
  // Compute sigmaTilde
  double* sigmaTilde = mat4;
  double* sigmaTemp = mat1; // Reuse mat1 as H is nolonger required
  // Compute F \widetilde{S}
  dgemm_(chn, chn, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
      STilde, &ndim, &zero, sigmaTemp, &ndim);
  double Jm53 = pow(Jm13, 5);
  // Compute F \widetilde{S} F^T
  // sigmaTilde = Ftilde Stilde Ftilde^T/J
  dgemm_(chn, chy, &ndim, &ndim, &ndim, &Jm53, sigmaTemp, &ndim,
      F_element_gp, &ndim, &zero, sigmaTilde, &ndim);

  // Compute deviatoric part of sigmaTilde and add the isochoric part to obtain
  // sigma
  const double traceSigmaTildeby3 = (sigmaTilde[0] + sigmaTilde[4] + sigmaTilde[8])/3.0;
  FILE_LOG_SINGLE(WARNING, "Trace of sigma Tilde = %15.9f", 3.0*traceSigmaTildeby3);
  const double volum = K*(J-1.0); 
  const double diagDiff = volum-traceSigmaTildeby3;
  // Compute sigma
  double* sigma_e = sigmaTilde;
  sigma_e[0] += diagDiff;
  sigma_e[4] += diagDiff;
  sigma_e[8] += diagDiff;

  /*
  //S = J*Finv*Sigma*FinvT
  double *STemp = mat3;//reuse mat3

  dgemm_(chn, chn, &ndim, &ndim, &ndim, &one, fInv, &ndim,
      sigma_e, &ndim, &zero, STemp, &ndim);
  dgemm_(chn, chy, &ndim, &ndim, &ndim, &J, STemp, &ndim,
      fInv, &ndim, &zero, S, &ndim);
  // FILE_LOGMatrix(WARNING, S, ndim, ndim, "Before viscoelaticty\n");
  double *Cmat = mat3;

  dgemm_(chy, chn, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
      F_element_gp, &ndim, &zero, Cmat, &ndim);
  // FILE_LOGMatrix(WARNING, Cmat, ndim, ndim, "C mat\n");
  // compute the double dot product C:S
  double SddC = 0.0;
  for (int i = 0; i < ndim2; ++i) {
    SddC += Cmat[i]*S[i];
  }
  // Compute isochoric part : Sic = 1/3(C:S)inv(C)
  SddC = SddC/3.0;
  // Compute Inverse of C, inv(C) = inv(F)*inv(F)^T
  double* Sic = mat3; // Reuse mat3 for storing isochoric part of S

  dgemm_(chn, chy, &ndim, &ndim, &ndim, &SddC, fInv, &ndim,
      fInv, &ndim, &zero, Sic, &ndim);
  // FILE_LOGMatrix(WARNING, Sic, ndim, ndim, "Sic mat\n");
  double *Sdev = mat2; // Reuse mat2 to store deviatoric part of S
  for (int i = 0; i < ndim2; ++i) {
    Sdev[i] = S[i]-Sic[i];
  }
  // FILE_LOGMatrix(WARNING, Sdev, ndim, ndim, "Sdev mat\n");
  // Access histroy dependence
  double *elemS = S0n[e];
  double *S0nLocal = &(elemS[gp*ndim2]);
  double *elemH = Hn[e];
  double *HnLocal = &(elemH[gp*ndim2*nPronyLocal]);
  for (unsigned int i = 0; i < nPronyLocal; ++i) {
    double *HnI = &(HnLocal[i*ndim2]);
    // Compute C1i = exp(-dt/t_i) and C2i = g_i (1-exp(-dt/t_i))/(dt/t_i)
    const double rt = dt/ti[i];
    const double c1i = exp(-rt);
    const double c2i = gi[i]*(1-c1i)/rt;
    // Update H_j^{n+1} = c1j*H_j^n + c2j*[Dev S_0^{n+1} - S_0^n]
    // Update S^{n+1} = S_0^{n+1} + \sum_j H_j^{n+1}
    for (int j = 0; j < ndim2; ++j) {
      HnI[j] = c1i*HnI[j]+c2i*(Sdev[j]-S0nLocal[j]);
      Sdev[j] = Sdev[j] + HnI[j];
    }
  }
  double *tau = mat4;
  dgemm_(chn, chn, &ndim, &ndim, &ndim, &one, F_element_gp, &ndim,
      Sdev, &ndim, &zero, STemp, &ndim);
  dgemm_(chn, chy, &ndim, &ndim, &ndim, &one, STemp, &ndim,
      F_element_gp, &ndim, &zero, tau, &ndim);
  tau[0] = tau[0] + J*K*(J-1);
  tau[4] = tau[4] + J*K*(J-1);
  tau[8] = tau[8] + J*K*(J-1);

  for(int i = 0; i< ndim2; ++i) {
    sigma_e[i] = tau[i]*invJ;
  }
  */

  // FILE_LOGMatrix(WARNING, S, ndim, ndim, "S Visco mat\n");

  // TODO (NEED UPDATE from PK2 to Cauchy)
  // Get location of array to store Cauchy values
  double * sigma_nLocal = &(sigma_n[sigmaptr[e]+6*gp]);
  // 6 values saved per gauss point for 3d
  // in voigt notation, sigma11
  sigma_nLocal[0] = sigma_e[0];
  // in voigt notation, sigma22
  sigma_nLocal[1] = sigma_e[4];
  // in voigt notation, sigma33
  sigma_nLocal[2] = sigma_e[8];
  // in voigt notation, sigma23
  sigma_nLocal[3] = sigma_e[7];
  // in voigt notation, sigma13
  sigma_nLocal[4] = sigma_e[6];
  // in voigt notation, sigma12
  sigma_nLocal[5] = sigma_e[3];

  // Update deviatoric part of S_0 for next time step
  // for (int i = 0; i < ndim2; ++i) {
  //   S0nLocal[i] = Sdev[i];
  // }
  free(gi);
  free(ti);

#ifdef DEBUG
  FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Element ID = %d, gp = %d", e, gp);
  for(int i=0;i<6;i++){
    int indexD = sigmaptr[e]+6*gp+i;
    FILE_LOG_SINGLE(DEBUGLOGIGNORE, "Sigma[%d] = %10.6e", indexD, sigma_n[indexD]);
  }
#endif //DEBUG
}
