#include"common.h"
#include"defs.h"
#include"externGlobal.h"
#include"util.h"
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
