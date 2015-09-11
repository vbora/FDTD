#include"common.h"
#include"defs.h"
#include"externGlobal.h"
#include"materials.h"
#include"dispersiveMat.h"
#include"util.h"
// Define all the materials here
void defMaterial(){

	//For each material define eps, mu, sigma
	int m, matIndex=0;
	int dm, dispMatIndex = 0;
	int lorentzIndex=0;
	int p;
	// material 0, air
	m = matIndex++;
	EPSR(m) = 1;	// constitutive parameters 
	MUR(m) = 1;	// constitutive parameters
	SIGMA(m) = 0;	// constitutive parameters
	JSX(m) = 0;
	JSY(m) = 0;
	JSZ(m) = 0;
	MSX(m) = 0;
	MSY(m) = 0;
	MSZ(m) = 0;
	setOrdinaryMaterial(m);
	printf("m %d, CA %e, CB, %e, DA %e, DB %e\n",m, CA(m), CB(m), DA(m), DB(m));

	// material 1, source
	m=matIndex++;
	EPSR(m) = 1;	// constitutive parameters 
	MUR(m) = 1;	// constitutive parameters
	SIGMA(m) = 0;	// constitutive parameters
	JSX(m) = 0;
	JSY(m) = 0;
	JSZ(m) = 0;
	MSX(m) = 0;
	MSY(m) = 0;
	MSZ(m) = 0;
	setOrdinaryMaterial(m);
	
	// make dielectric
	
	// material 2, dielectric 
	m = matIndex++;
	EPSR(m) = 5;	// constitutive parameters 
	MUR(m) = 1;	// constitutive parameters
	SIGMA(m) = 0;	// constitutive parameters
	JSX(m) = 0;
	JSY(m) = 0;
	JSZ(m) = 0;
	MSX(m) = 0;
	MSY(m) = 0;
	MSZ(m) = 0;
	setOrdinaryMaterial(m);
	printf("m %d, CA %e, CB, %e, DA %e, DB %e\n",m, CA(m), CB(m), DA(m), DB(m));

	if (PML_LAYERS > 0)
		defPMLMaterial();

	// material 3, 
	m = matIndex++;
	NPOLES(m) = 2; // Two pole
	// Only give the value of lorentzIndex
	EPSR(m) = 3.7221; // eps-infinity 
	MUR(m) = 1;	// constitutive parameters
	SIGMA(m) = 0;	// constitutive parameters
	JSX(m) = 0;
	JSY(m) = 0;
	JSZ(m) = 0;
	MSX(m) = 0;
	MSY(m) = 0;
	MSZ(m) = 0;
	LORN_INDEX(m) = lorentzIndex;
	lorentzIndex += NPOLES(m);
	// Now specify the  parameter for each pole
	p=LORN_INDEX(m);

	DELTAEPS(p) = 1.0524;
	OMEGAP(p)   = 0.0241;
	DELTAP(p)   = 0.0032;

	p++;

	DELTAEPS(p) = 0.0048e6;
	OMEGAP(p)   = 1e-6;
	DELTAP(p)   = 8.142e-5;
			
	defDispersiveConstants(m);//Calculate the constantsi

	totalNumberOfMaterials = matIndex;			
} // end defMaterial 



// Define geometry here
void defGeometry(){
	int m;
	int ix,iy,iz,gx,gy,gz;
	int center1[3], radius1; 
	int center2[3], radius2;

	//Ca,Cb is a function of sigma, deltaT, epsilon,
	//Da. Db is a function of sigma, deltaT, mu  

	// the xLen, yLen also include the extra lines for boundary conditions
	for (ix=0; ix < xLen + 2; ix++)
	for (iy=0; iy < yLen + 2; iy++)
	for (iz=0; iz < zLen + 2; iz++){
		MATERIALINDEX(ix,iy,iz)=3;
	}// end all for
	
	// set source
	m = 1 ;

	//NOTE these are global coordinates
	// Source travelling in Y direction
	/*gy = 1;
	for (gx=1;gx<XDIM;gx++)
		for (gz=1;gz<ZDIM;gz++)
			if (isLocal(gx,gy,gz,&ix,&iy,&iz))		
 		       		MATERIALINDEX(ix,iy,iz)=m;
				
	*/

	// Source travelling in X direction
	//XDIM, YDIM are the total simulation dimensions
	gx = PML_LAYERS+15;
	for (gy=0;gy<YDIM;gy++)
		for (gz=0;gz<ZDIM;gz++)
			if (isLocal(gx,gy,gz,&ix,&iy,&iz)){
				//printf("source processor sid # %d  gx %d gy %d gz %d ix %d iy %d  iz  %d\n",sid,gx,gy,gz,ix,iy,iz);	
 		       		MATERIALINDEX(ix, iy,iz)=m;
			}
/*	
	m = 3;

	for (ix=0; ix < xLen + 2; ix++)
	for (iy=0; iy < yLen + 2; iy++)
	for (iz=0; iz < zLen + 2; iz++){
		globalEquivalent(ix, iy, iz, &gx, &gy, &gz);
		if ( gx>XDIM/2)
			MATERIALINDEX(ix,iy,iz)=3;
	}// end all for
*/
/*
	// Make sphere
	m = 3;
	radius1 = radius2 = 12.5;	

	center1[X] = XDIM/2 - 1 ; 
	center1[Y] = YDIM/2 - 1 - 25;		
	center1[Z] = ZDIM/2 - 1;		
	makeSphere(m,center1,radius1);
 
	center2[X] = XDIM/2 - 1 ; 
	center2[Y] = YDIM/2 - 1 + 25;		
	center2[Z] = ZDIM/2 - 1 ; 
	makeSphere(m,center2,radius2);
*/
//defining PML layers
// for now, just putting in a huge loop, but it if really affects the performance, then we have to do a piecewise solution

}// End defGeometry

// Source parameters here

void defSources(){
	int m;
	int i=0; // source number
	// source 1
	/*
	m = 1;
	listOfSources[0] = m;
	SOURCETYPE(m) = PULSE ;
	DELAY(m) =15 ;
	POLARIZATION(m)= HZFIELD;
	PEAKVALUE(m) = 100; 
	LOCATION(m) =1;
	*/
	// 1nm --> 5 points
	m = 1;
	listOfSources[i++] = m;
	SOURCETYPE(m) = SINE; 
	FREQUENCY(m)  = 0.0018762; // 1/533; 
	//DELTA_FREQUENCY_M not required 
	PHASE(m) = 0; //in terms of pi 
	DELAY(m) = 0; // time steps 
	POLARIZATION(m) = HZFIELD; 
	PEAKVALUE(m) = 1; 
	LOCATION(m) = 1; 

	listOfSources[i] = -1; // terminating condition
}
