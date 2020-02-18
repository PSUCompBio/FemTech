#include "FemTech.h"

double volumeHexahedron(double *Lcoordinates) {
  double q0[3], q1[3], q2[3], q3[3], q4[3], q5[3], q6[3], q7[3];
  for (int i = 0; i < 3; ++i) {
    q0[i] = Lcoordinates[3*0+i]-Lcoordinates[3*1+i]+Lcoordinates[3*2+i]-\
            Lcoordinates[3*3+i]+Lcoordinates[3*4+i]-Lcoordinates[3*5+i]+\
            Lcoordinates[3*6+i]-Lcoordinates[3*7+i];
    q1[i] = Lcoordinates[3*0+i]-Lcoordinates[3*1+i]-Lcoordinates[3*2+i]+\
            Lcoordinates[3*3+i]-Lcoordinates[3*4+i]+Lcoordinates[3*5+i]+\
            Lcoordinates[3*6+i]-Lcoordinates[3*7+i];
    q2[i] = -Lcoordinates[3*0+i]+Lcoordinates[3*1+i]+Lcoordinates[3*2+i]-\
            Lcoordinates[3*3+i]-Lcoordinates[3*4+i]+Lcoordinates[3*5+i]+\
            Lcoordinates[3*6+i]-Lcoordinates[3*7+i];
    q3[i] = Lcoordinates[3*0+i]+Lcoordinates[3*1+i]-Lcoordinates[3*2+i]-\
            Lcoordinates[3*3+i]-Lcoordinates[3*4+i]-Lcoordinates[3*5+i]+\
            Lcoordinates[3*6+i]+Lcoordinates[3*7+i];
    q4[i] = -Lcoordinates[3*0+i]-Lcoordinates[3*1+i]+Lcoordinates[3*2+i]+\
            Lcoordinates[3*3+i]-Lcoordinates[3*4+i]-Lcoordinates[3*5+i]+\
            Lcoordinates[3*6+i]+Lcoordinates[3*7+i];
    q5[i] = -Lcoordinates[3*0+i]-Lcoordinates[3*1+i]-Lcoordinates[3*2+i]-\
            Lcoordinates[3*3+i]+Lcoordinates[3*4+i]+Lcoordinates[3*5+i]+\
            Lcoordinates[3*6+i]+Lcoordinates[3*7+i];
  }
  return (tripleProduct(q0, q4, q3) + tripleProduct(q2, q0, q1) +
      tripleProduct(q1, q3, q5))/192.0 + tripleProduct(q2, q4, q5)/64.0;
}

double areaHexahedronFace(double *Lcoordinates, const int * const index) {
  // Points are assumed in order
  // Store the points from index
  double *p1 = &(Lcoordinates[3*index[0]]);
  double *p2 = &(Lcoordinates[3*index[1]]);
  double *p3 = &(Lcoordinates[3*index[2]]);
  double *p4 = &(Lcoordinates[3*index[3]]);

  double tol = 1e-6;
  double centerD[3];
  double c1[3], c2[3];
  for (int i = 0; i < 3; ++i) {
    centerD[i] = 0.25*(p1[i]-p2[i]+p3[i]-p4[i]);
    c1[i] = 0.25*(-p1[i]+p2[i]+p3[i]-p4[i]);
    c2[i] = 0.25*(-p1[i]-p2[i]+p3[i]+p4[i]);
  }
  // Check if its a parallelogram
  if ((centerD[0] < tol) && (centerD[1] < tol) && (centerD[2] < tol)) {
    return 4.0*normOfCrossProduct(c1, c2);
  }
  // Otherwise use 2x2 Quadrature
  double t = sqrt(3.0)/3.0;
  const double q[2] = {-t, t};
  double area = 0.0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      double v1[3], v2[3];
      for (int k = 0; k < 3; ++k) {
        v1[k] = q[j]*centerD[k]+c1[k]; 
        v2[k] = q[i]*centerD[k]+c2[k]; 
      }
      area += normOfCrossProduct(v1, v2);
    }
  }
  return area;
}

double volumeTetrahedron(double *Lcoordinates) {
  double a[3], b[3], c[3];
  double volume;
  for (int i = 0; i < 3; ++i) {
    a[i] = Lcoordinates[i]-Lcoordinates[9+i];
    b[i] = Lcoordinates[3+i]-Lcoordinates[9+i];
    c[i] = Lcoordinates[6+i]-Lcoordinates[9+i];
  }
  volume = fabs(tripleProduct(a, b, c))/6.0;
  return volume;
}
