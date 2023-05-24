#include "FemTech.h"

void FindNaturalCoord(){
   // Purpose : Compute chi, eta, iota for fiber nodes
  int hostel;
  double a[3], xs[9], f[3], J_Inv[9], detJembed, error[3], maxerrortemp, incr[3];
  double tol = 1e-6;
  double N1, N2, N3, N4, N5, N6, N7, N8;
  double DN1Dchi, DN2Dchi, DN3Dchi, DN4Dchi, DN5Dchi, DN6Dchi, DN7Dchi, DN8Dchi, DN1Deta, DN2Deta, DN3Deta, DN4Deta, DN5Deta, DN6Deta, DN7Deta, DN8Deta, DN1Diota, DN2Diota, DN3Diota, DN4Diota, DN5Diota, DN6Diota, DN7Diota, DN8Diota;
  int n_host = nelements - nembedel;
  for(int i = 0; i<nNodes; i++){
    double maxerror = 1;
     if(nodeconstrain[i]!=-1){
	int hostelglobal = nodeconstrain[i];
	for(int j1 = 0; j1<nelements; j1++){
		if(hostelglobal+1==global_eid[j1]){
			hostel = j1;
			a[0] = coordinates[i*ndim + 0];
			a[1] = coordinates[i*ndim + 1];
			a[2] = coordinates[i*ndim + 2]; 
			double chi, eta, iota, chinew, etanew, iotanew;
			chi =  0.1;
			eta =  0.1;
			iota = 0.1;
			N1 = ((1 - chi)*(1 - eta)*(1 - iota)) / 8;
			N2 = ((1 + chi)*(1 - eta)*(1 - iota)) / 8;
			N3 = ((1 + chi)*(1 + eta)*(1 - iota)) / 8;
			N4 = ((1 - chi)*(1 + eta)*(1 - iota)) / 8;
			N5 = ((1 - chi)*(1 - eta)*(1 + iota)) / 8;
			N6 = ((1 + chi)*(1 - eta)*(1 + iota)) / 8;
			N7 = ((1 + chi)*(1 + eta)*(1 + iota)) / 8;
			N8 = ((1 - chi)*(1 + eta)*(1 + iota)) / 8;
			while(maxerror>tol){
				// The first derivatives with respect to chi
				DN1Dchi = -((eta - 1)*(iota - 1)) / 8;
				DN2Dchi =  ((eta - 1)*(iota - 1)) / 8;
				DN3Dchi = -((eta + 1)*(iota - 1)) / 8;
				DN4Dchi =  ((eta + 1)*(iota - 1)) / 8;
				DN5Dchi =  ((eta - 1)*(iota + 1)) / 8;
				DN6Dchi = -((eta - 1)*(iota + 1)) / 8;
				DN7Dchi =  ((eta + 1)*(iota + 1)) / 8;
				DN8Dchi = -((eta + 1)*(iota + 1)) / 8;
				// with respect to eta
				DN1Deta = -((chi - 1)*(iota - 1)) / 8;
				DN2Deta =  ((chi + 1)*(iota - 1)) / 8;
				DN3Deta = -((chi + 1)*(iota - 1)) / 8;
				DN4Deta =  ((chi - 1)*(iota - 1)) / 8;
				DN5Deta =  ((chi - 1)*(iota + 1)) / 8;
				DN6Deta = -((chi + 1)*(iota + 1)) / 8;
				DN7Deta =  ((chi + 1)*(iota + 1)) / 8;
				DN8Deta = -((chi - 1)*(iota + 1)) / 8;
				// with respect to iota
				DN1Diota = -((chi - 1)*(eta - 1)) / 8;
				DN2Diota =  ((chi + 1)*(eta - 1)) / 8;
				DN3Diota = -((chi + 1)*(eta + 1)) / 8;
				DN4Diota =  ((chi - 1)*(eta + 1)) / 8;
				DN5Diota =  ((chi - 1)*(eta - 1)) / 8;
				DN6Diota = -((chi + 1)*(eta - 1)) / 8;
				DN7Diota =  ((chi + 1)*(eta + 1)) / 8;
				DN8Diota = -((chi - 1)*(eta + 1)) / 8;
				//compute the Jacobian inverse
				for (int j = 0; j<ndim; j++){		
					xs[j+0*ndim] = DN1Dchi*coordinates[ndim*connectivity[eptr[hostel]+0] + j] + DN2Dchi*coordinates[ndim*connectivity[eptr[hostel]+1] + j] + DN3Dchi*coordinates[ndim*connectivity[eptr[hostel]+2] + j]
						+ DN4Dchi*coordinates[ndim*connectivity[eptr[hostel]+3] + j] + DN5Dchi*coordinates[ndim*connectivity[eptr[hostel]+4] + j] + DN6Dchi*coordinates[ndim*connectivity[eptr[hostel]+5] + j]
						+ DN7Dchi*coordinates[ndim*connectivity[eptr[hostel]+6] + j] + DN8Dchi*coordinates[ndim*connectivity[eptr[hostel]+7] + j];

					xs[j+1*ndim] = DN1Deta*coordinates[ndim*connectivity[eptr[hostel]+0] + j] + DN2Deta*coordinates[ndim*connectivity[eptr[hostel]+1] + j] + DN3Deta*coordinates[ndim*connectivity[eptr[hostel]+2] + j]
						+ DN4Deta*coordinates[ndim*connectivity[eptr[hostel]+3] + j] + DN5Deta*coordinates[ndim*connectivity[eptr[hostel]+4] + j] + DN6Deta*coordinates[ndim*connectivity[eptr[hostel]+5] + j]
						+ DN7Deta*coordinates[ndim*connectivity[eptr[hostel]+6] + j] + DN8Deta*coordinates[ndim*connectivity[eptr[hostel]+7] + j];

					xs[j+2*ndim] = DN1Diota*coordinates[ndim*connectivity[eptr[hostel]+0] + j] + DN2Diota*coordinates[ndim*connectivity[eptr[hostel]+1] + j] + DN3Diota*coordinates[ndim*connectivity[eptr[hostel]+2] + j]
						+ DN4Diota*coordinates[ndim*connectivity[eptr[hostel]+3] + j] + DN5Diota*coordinates[ndim*connectivity[eptr[hostel]+4] + j] + DN6Diota*coordinates[ndim*connectivity[eptr[hostel]+5] + j]
						+ DN7Diota*coordinates[ndim*connectivity[eptr[hostel]+6] + j] + DN8Diota*coordinates[ndim*connectivity[eptr[hostel]+7] + j];
					}
			  	detJembed = inverse3x3Matrix(xs, J_Inv);
				//compute f = N*x-a
			  	for(int j=0; j<ndim; j++){
					  f[j] = N1*coordinates[ndim*connectivity[eptr[hostel]+0] + j] + N2*coordinates[ndim*connectivity[eptr[hostel]+1] + j] + N3*coordinates[ndim*connectivity[eptr[hostel]+2] + j] +
						 N4*coordinates[ndim*connectivity[eptr[hostel]+3] + j] + N5*coordinates[ndim*connectivity[eptr[hostel]+4] + j] + N6*coordinates[ndim*connectivity[eptr[hostel]+5] + j] + 
						 N7*coordinates[ndim*connectivity[eptr[hostel]+6] + j] + N8*coordinates[ndim*connectivity[eptr[hostel]+7] + j] - a[j];
				}
				//compute increment
				for(int j=0; j<ndim; j++){
					incr[j] = -(J_Inv[0*ndim+j]*f[0] + J_Inv[1*ndim+j]*f[1] + J_Inv[2*ndim+j]*f[2]);
				}
				chinew = incr[0] + chi;
				etanew = incr[1] + eta;
				iotanew = incr[2] + iota;
				chi = chinew;
				eta = etanew;
				iota = iotanew;
				N1 = ((1 - chi)*(1 - eta)*(1 - iota)) / 8;
				N2 = ((1 + chi)*(1 - eta)*(1 - iota)) / 8;
				N3 = ((1 + chi)*(1 + eta)*(1 - iota)) / 8;
				N4 = ((1 - chi)*(1 + eta)*(1 - iota)) / 8;
				N5 = ((1 - chi)*(1 - eta)*(1 + iota)) / 8;
				N6 = ((1 + chi)*(1 - eta)*(1 + iota)) / 8;
				N7 = ((1 + chi)*(1 + eta)*(1 + iota)) / 8;
				N8 = ((1 - chi)*(1 + eta)*(1 + iota)) / 8;
				for(int j=0; j<ndim; j++){
					error[j] = N1*coordinates[ndim*connectivity[eptr[hostel]+0] + j] + N2*coordinates[ndim*connectivity[eptr[hostel]+1] + j] + N3*coordinates[ndim*connectivity[eptr[hostel]+2] + j] +
						N4*coordinates[ndim*connectivity[eptr[hostel]+3] + j] + N5*coordinates[ndim*connectivity[eptr[hostel]+4] + j] + N6*coordinates[ndim*connectivity[eptr[hostel]+5] + j] + 
						N7*coordinates[ndim*connectivity[eptr[hostel]+6] + j] + N8*coordinates[ndim*connectivity[eptr[hostel]+7] + j] - a[j];
					}
				maxerrortemp = abs(error[0])>abs(error[1])?abs(error[0]):abs(error[1]);
				maxerror = maxerrortemp>abs(error[2])?maxerrortemp:abs(error[2]);
			  }
		  embedNC[i*ndim+0] = chi;
		  embedNC[i*ndim+1] = eta;
		  embedNC[i*ndim+2] = iota;
		  for(int j = 0; j<ndim; j++){
			if(abs(embedNC[i*ndim+j])<1e-12)
				embedNC[i*ndim+j] = 0.0;
			}
		  break;
		  }
	    }
       }
  }
  return;
}
