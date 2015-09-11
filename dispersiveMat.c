#include"common.h"
#include"defs.h"
#include "externGlobal.h"
#include "util.h"
#include "dispersiveMat.h"

// Make all the constants
void defDispersiveConstants(int m){

int i,numberOfPoles, index, p;
double sumGamma=0; 

index = LORN_INDEX(m); // Location of first pole
numberOfPoles = NPOLES(m); 

	for (i=0; i<numberOfPoles; i++){
		p = index+i;

		ALPHA(p) = (2 - OMEGAP(p)*OMEGAP(p)*deltaT*deltaT) / (1 + DELTAP(p)*deltaT);
		ZETA(p) = (DELTAP(p)*deltaT - 1 )  / (DELTAP(p)*deltaT + 1);
		GAMMA(p) = EPSNOT*DELTAEPS(p)*OMEGAP(p)*OMEGAP(p)*deltaT*deltaT / (DELTAP(p)*deltaT + 1);
	
		sumGamma +=GAMMA(p); 	
		//printf("processor # %d  p %d, ALPHA(p) %e ZETA(p) %e GAMMA(p) %e \n", sid, p, ALPHA(p), ZETA(p), GAMMA(p));
	} // end for

	DISPC1(m) = 0.5* sumGamma / (2*EPSNOT*EPSR(m) + 0.5*sumGamma + SIGMA(m)*deltaT);
	DISPC2(m) = (2*EPSNOT*EPSR(m) - SIGMA(m)*deltaT) / (2*EPSNOT*EPSR(m) + 0.5*sumGamma + SIGMA(m)*deltaT);
	DISPC3(m) = 0.5* deltaT / (2*EPSNOT*EPSR(m) + 0.5*sumGamma + SIGMA(m)*deltaT);

	//printf("processor # %d  DISPC1(m) %e,  DISPC2(m) %e,  DISPC3(m) %e \n", sid, DISPC1(m),  DISPC2(m),  DISPC3(m));

} //defDispersiveConstants

void allocateDispFields(){
	int dispMat[MAXMATERIALS], numDispMat=0;
	int count = 0;	
	int i,j,k,m,d;
	// Find out the number of dispersive materials	
	for (m=0; m<totalNumberOfMaterials; m++){
		if (NPOLES(m) != 0)
			dispMat[numDispMat++] = m;
	}

	// We populate the dispField array here
	for (i=1; i < xLen + 1; i++)
        for (j=1; j < yLen + 1; j++)
	for (k=1; k < zLen + 1; k++)
	for (d=0; d < numDispMat; d++){
		if (MATERIALINDEX(i,j,k) == dispMat[d]){
			DISPFIELDINDEX(i,j,k) = count; 
			count += NPOLES(dispMat[d]); 	
		}		
	}// end for

	dispField = (dispFieldStruct*)calloc(count, sizeof(dispFieldStruct));
	printf("sid # %d count in dispField %d numDispMat %d \n", sid, count, numDispMat);
	if (dispField == NULL){
		errorFunction("Malloc failed in dispField");
	}
} //end allocateDispFields	

void updateDispField(int i, int j, int k, int m, int n){
	
	int d,nPoles,p;	
	int lorentzIndex;
	double exN, eyN, ezN;
	double jdxN, jdyN, jdzN;
	double sumJdx, sumJdy, sumJdz;

	d = DISPFIELDINDEX(i,j,k);
	
	lorentzIndex = LORN_INDEX(m);
	
	nPoles = NPOLES(m);
	exN = EX(i,j,k);
	eyN = EY(i,j,k);
	ezN = EZ(i,j,k);

	// for calculating E	
	for (p=0; p<nPoles; p++){
		sumJdx = (1+ALPHA(lorentzIndex+p)) * JDX(d+p) + ZETA(lorentzIndex+p) * JDXMINUS1(d+p);
		sumJdy = (1+ALPHA(lorentzIndex+p)) * JDY(d+p) + ZETA(lorentzIndex+p) * JDYMINUS1(d+p);
		sumJdz = (1+ALPHA(lorentzIndex+p)) * JDZ(d+p) + ZETA(lorentzIndex+p) * JDZMINUS1(d+p);
	}// end for
		
		
	EX(i,j,k)= DISPC1(m) * EXMINUS1(d) + DISPC2(m) *  EX(i,j,k) + DISPC3(m) * ( HZ(i,j,k) - HZ(i,j-1,k) + HY(i,j,k-1) - HY(i,j,k) - 0.5 * sumJdx);
	EY(i,j,k)= DISPC1(m) * EYMINUS1(d) + DISPC2(m) *  EY(i,j,k) + DISPC3(m) * ( HX(i,j,k) - HX(i,j,k-1) + HZ(i-1,j,k) - HZ(i,j,k) - 0.5 * sumJdy);
	EZ(i,j,k)= DISPC1(m) * EZMINUS1(d) + DISPC2(m) *  EZ(i,j,k) + DISPC3(m) * ( HY(i,j,k) - HY(i-1,j,k) + HX(i,j-1,k) - HX(i,j,k) - 0.5 * sumJdz);

	
	// After E updating the currents	
	for (p=0; p<nPoles; p++){
		jdxN =  JDX(d+p);
		jdyN =  JDY(d+p);
		jdzN =  JDZ(d+p);

		JDX(d+p) = ALPHA(lorentzIndex+p) * JDX(d+p) + ZETA(lorentzIndex+p) * JDXMINUS1(d+p) + GAMMA(lorentzIndex+p) * ((EX(i,j,k)-EXMINUS1(d)) / (2* deltaT)) ;

		JDY(d+p) = ALPHA(lorentzIndex+p) * JDY(d+p) + ZETA(lorentzIndex+p) * JDYMINUS1(d+p) + GAMMA(lorentzIndex+p) * ((EY(i,j,k)-EYMINUS1(d)) / (2* deltaT)) ;

		JDZ(d+p) = ALPHA(lorentzIndex+p) * JDZ(d+p) + ZETA(lorentzIndex+p) * JDZMINUS1(d+p) + GAMMA(lorentzIndex+p) * ((EX(i,j,k)-EXMINUS1(d)) / (2* deltaT)) ;

		// Now we have to update the values in JDXMINUS1
		JDXMINUS1(d+p) = jdxN;
		JDYMINUS1(d+p) = jdyN;
		JDZMINUS1(d+p) = jdzN;

	}// end for
		
	// And update EXMINUS1
	EXMINUS1(d) = exN;
	EYMINUS1(d) = eyN;
	EZMINUS1(d) = ezN;
	
	//printf("\n processor sid %d TIMESTEP %d EX(i,j,k) %e, exN %e, EY(i,j,k) %e, eyN %e, EZ(i,j,k) %e, ezN %e ",sid, n, EX(i,j,k), exN, EY(i,j,k), eyN, EZ(i,j,k), ezN );	
} // end function
