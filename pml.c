#include"common.h"
#include"defs.h"
#include"externGlobal.h"
#include"util.h"
#include <stdlib.h>

void allocatePMLFields(){
	int i,j,k;
	int gx,gy,gz;
	long count;
	count=0;
	for (i=1;i<xLen+1; i++)
	for (j=1;j<yLen+1; j++)
	for (k=1;k<zLen+1; k++){
		globalEquivalent(i,j,k, &gx,&gy,&gz);
		if ((gx<PML_LAYERS)||(XDIM-(gx+1) < PML_LAYERS) || (gy<PML_LAYERS)||(YDIM-(gy+1) < PML_LAYERS) || (gz<PML_LAYERS)||(ZDIM-(gz+1) < PML_LAYERS))	
			PMLINDEX(i,j,k) = count++;	
				
	} // end for
	printf("sid %d count %ld\n", sid, count);
	pmlFields = (pmlFieldStruct *)calloc(count, sizeof(pmlFieldStruct));
	if (pmlFields == NULL){
		errorFunction("Malloc failed in pml fields");
        }
}
	

void defPMLMaterial(){
        int m,i, index;
        double g,delta;
        double sigma, sigmaMax, kappa;

        sigmaMax = -(EXP_M + 1)* log(REFL0)/(2*ETTA*PML_LAYERS);

        for (i=0; i<PML_LAYERS; i++){
                sigma = pow(((double)i+1 / (PML_LAYERS)),EXP_M)* sigmaMax;
                kappa = 1 + (KAPPA_MAX -1 )*pow( ((double)i+1/PML_LAYERS), EXP_M);    

                cn[PML_LAYERS-1-i].c1 = (2*EPSNOT*kappa - sigma*deltaT) / (2*EPSNOT*kappa + sigma*deltaT);
                cn[PML_LAYERS-1-i].c2 = (2*EPSNOT*deltaT) / (2*EPSNOT*kappa + sigma*deltaT);  
                cn[PML_LAYERS-1-i].c3 = cn[PML_LAYERS-i-1].c1; 
                cn[PML_LAYERS-1-i].c4 = 1 / (2*EPSNOT*kappa + sigma*deltaT);
                cn[PML_LAYERS-1-i].c5 = (2*EPSNOT*kappa + sigma*deltaT);
                cn[PML_LAYERS-1-i].c6 = (2*EPSNOT*kappa - sigma*deltaT);
		if (sid == 0)
			printf("i %d sigma %e sigmaMax %e kappa %e c1 %e, c2 %e, c3 %e, c4 %e, c5 %e, c6 %e \n",i,sigma, sigmaMax, kappa,  cn[PML_LAYERS-1-i].c1,cn[PML_LAYERS-1-i].c2, cn[PML_LAYERS-1-i].c3, cn[PML_LAYERS-1-i].c4, cn[PML_LAYERS-1-i].c5, cn[PML_LAYERS-1-i].c6 );

                cOnes[i].c1 = 1;
                cOnes[i].c2 = deltaT;
                cOnes[i].c3 = 1;
                cOnes[i].c4 = 1/(2*EPSNOT);
                cOnes[i].c5 = (2*EPSNOT);
                cOnes[i].c6 = (2*EPSNOT);
//		if (sid == 0)
//			printf("cOnes i %d c1 %e, c2 %e, c3 %e, c4 %e, c5 %e, c6 %e \n",i,cOnes[i].c1,cOnes[i].c2, cOnes[i].c3, cOnes[i].c4, cOnes[i].c5, cOnes[i].c6 );
	}
} //end defPMLMaterial
/*
void defPMLMaterial(){
        int m,i, index;
        double g,delta;
        double sigma, sigmaMax, kappa;

        sigmaMax = -(EXP_M + 1)* log(REFL0)/(2*ETTA*PML_LAYERS);

        for (i=0; i<PML_LAYERS; i++){
                sigma = pow(((double)i / (PML_LAYERS-1)),EXP_M)* sigmaMax;
                kappa = 1 + (KAPPA_MAX -1 )*powf( (i/PML_LAYERS), EXP_M);    
                cInv[i].c1 = (2*EPSNOT*kappa - sigma*deltaT) / (2*EPSNOT*kappa + sigma*deltaT);
                cInv[i].c2 = (2*EPSNOT*deltaT) / (2*EPSNOT*kappa + sigma*deltaT);
                cInv[i].c3 = cInv[i].c1;
                cInv[i].c4 = 1 / (2*EPSNOT*kappa + sigma*deltaT);
                cInv[i].c5 = (2*EPSNOT*kappa + sigma*deltaT);
                cInv[i].c6 = (2*EPSNOT*kappa - sigma*deltaT);

                cn[PML_LAYERS-i-1].c1 =  cInv[i].c1;
                cn[PML_LAYERS-i-1].c2 =  cInv[i].c2;
                cn[PML_LAYERS-i-1].c3 =  cInv[i].c3;
                cn[PML_LAYERS-i-1].c4 =  cInv[i].c4;
                cn[PML_LAYERS-i-1].c5 =  cInv[i].c5;
                cn[PML_LAYERS-i-1].c6 =  cInv[i].c6;

		printf("i %d sigma %e sigmaMax %e c1 %e, c2 %e, c3 %e, c4 %e, c5 %e, c6 %e \n",i,sigma, sigmaMax, cInv[i].c1,cInv[i].c2, cInv[i].c3, cInv[i].c4, cInv[i].c5, cInv[i].c6 );

                cOnes[i].c1 = 1;
                cOnes[i].c2 = deltaT;
                cOnes[i].c3 = 1;
                cOnes[i].c4 = 1/(2*EPSNOT);
                cOnes[i].c5 = (2*EPSNOT);
                cOnes[i].c6 = (2*EPSNOT);
		printf("cOnes i %d c1 %e, c2 %e, c3 %e, c4 %e, c5 %e, c6 %e \n",i,cOnes[i].c1,cOnes[i].c2, cOnes[i].c3, cOnes[i].c4, cOnes[i].c5, cOnes[i].c6 );
	}
} //end defPMLMaterial
*/
