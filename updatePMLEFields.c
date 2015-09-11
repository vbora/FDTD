#include"common.h"
#include"defs.h"
#include"externGlobal.h"
#include"util.h"
#include"updateFields.h" 
#include<stdlib.h>
#include"updatePMLFields.h"

void updatePMLEFieldUtilAll(int xStart, int xEnd, int yStart, int yEnd, int zStart, int zEnd, int xBound, int yBound, int zBound, pmlConstStruct cx[], pmlConstStruct cy[], pmlConstStruct cz[], int n){
		
	int i, j, k; // These are the indexes
	int x,y,z;// There are indexes for the pmlConstStruct indexs
	long p;
	int m;
	double dxStore, dyStore, dzStore;
	double eps;
// Either this, or make cOnes a bigger array with number of elements which are the max of xLen, yLen, zLen
	for (i=xStart;i<xEnd+1;i++){
		x = abs(xBound - i)%PML_LAYERS;	
	for (j=yStart;j<yEnd+1;j++){
		y = abs(yBound - j)%PML_LAYERS;		
	for (k=zStart;k<zEnd+1;k++){
		z = abs(zBound - k)%PML_LAYERS;
			
		//printf("(i,j,k) (%d, %d, %d), (x,y,z) (%d,%d,%d)\n", i,j,k,x,y,z);
		p = PMLINDEX(i,j,k);

		dxStore = DX(p);					
		dyStore = DY(p);					
		dzStore = DZ(p);			
	
		m = MATERIALINDEX(i,j,k);
		eps = EPSR(m)*EPSNOT;
				
		DX(p)     = cy[y].c1*DX(p) + cy[y].c2*(HZ(i,j,k) - HZ(i,j-1,k) - HY(i,j,k) + HY(i,j,k-1) - JSX(m));
		EX(i,j,k) = cz[z].c3*EX(i,j,k) + cz[z].c4*(cx[x].c5*DX(p)-cx[x].c6*dxStore)/eps; 
		 
		DY(p)     = cz[z].c1*DY(p) + cz[z].c2*(HX(i,j,k) - HX(i,j,k-1) + HZ(i-1,j,k) - HZ(i,j,k) - JSY(m));
		EY(i,j,k) = cx[x].c3*EY(i,j,k) + cx[x].c4*(cy[y].c5*DY(p)-cy[y].c6*dyStore)/eps; 
		
		DZ(p)     = cx[x].c1*DZ(p) + cx[x].c2*(HY(i,j,k) - HY(i-1,j,k) - HX(i,j,k) + HX(i,j-1,k) - JSZ(m));
		EZ(i,j,k) = cy[y].c3*EZ(i,j,k) + cy[y].c4*(cz[z].c5*DZ(p)-cz[z].c6*dzStore)/eps; 
	}}} // end for

}// end updatePMLEFieldUtilAll

void updatePMLEFieldUtilEx(int xStart, int xEnd, int yStart, int yEnd, int zStart, int zEnd, int xBound, int yBound, int zBound, pmlConstStruct cx[], pmlConstStruct cy[], pmlConstStruct cz[], int n){
		
	int i, j, k; // These are the indexes
	int x,y,z;// There are indexes for the pmlConstStruct indexs
	long p;
	int m;
	double dxStore;
	double eps;
	for (i=xStart;i<xEnd+1;i++){
		x = abs(xBound - i)%PML_LAYERS;	
	for (j=yStart;j<yEnd+1;j++){
		y = abs(yBound - j)%PML_LAYERS;	
	for (k=zStart;k<zEnd+1;k++){
		z = abs(zBound - k)%PML_LAYERS;
		
		p = PMLINDEX(i,j,k);
		dxStore = DX(p);					
	
		m = MATERIALINDEX(i,j,k);
		eps = EPSR(m)*EPSNOT;
		DX(p)     = cy[y].c1*DX(p) + cy[y].c2*(HZ(i,j,k) - HZ(i,j-1,k) - HY(i,j,k) + HY(i,j,k-1) - JSX(m));
		EX(i,j,k) = cz[z].c3*EX(i,j,k) + cz[z].c4*(cx[x].c5*DX(p)-cx[x].c6*dxStore)/eps; 
		 
	}}} // end for

}// end updatePMLEFieldUtilEx

void updatePMLEFieldUtilEy(int xStart, int xEnd, int yStart, int yEnd, int zStart, int zEnd, int xBound, int yBound, int zBound, pmlConstStruct cx[], pmlConstStruct cy[], pmlConstStruct cz[], int n){
		
	int i, j, k; // These are the indexes
	int x,y,z;// There are indexes for the pmlConstStruct indexs
	long p;
	int m;
	double dyStore;
	double eps;
	for (i=xStart;i<xEnd+1;i++){
		x = abs(xBound - i)%PML_LAYERS;	
	for (j=yStart;j<yEnd+1;j++){
		y = abs(yBound - j)%PML_LAYERS;	
	for (k=zStart;k<zEnd+1;k++){
		z = abs(zBound - k)%PML_LAYERS;
			
		p = PMLINDEX(i,j,k);
		dyStore = DY(p);					
	
		m = MATERIALINDEX(i,j,k);
		eps = EPSR(m)*EPSNOT;
				
		DY(p)     = cz[z].c1*DY(p) + cz[z].c2*(HX(i,j,k) - HX(i,j,k-1) + HZ(i-1,j,k) - HZ(i,j,k) - JSY(m));
		EY(i,j,k) = cx[x].c3*EY(i,j,k) + cx[x].c4*(cy[y].c5*DY(p)-cy[y].c6*dyStore)/eps; 

	}}} // end for

}// end updatePMLEFieldUtilEy

void updatePMLEFieldUtilEz(int xStart, int xEnd, int yStart, int yEnd, int zStart, int zEnd, int xBound, int yBound, int zBound, pmlConstStruct cx[], pmlConstStruct cy[], pmlConstStruct cz[], int n){
		
	int i, j, k; // These are the indexes
	int x,y,z;// There are indexes for the pmlConstStruct indexs
	long p;
	int m;
	double dzStore;
	double eps;
	for (i=xStart;i<xEnd+1;i++){
		x = abs(xBound - i)%PML_LAYERS;	
	for (j=yStart;j<yEnd+1;j++){
		y = abs(yBound - j)%PML_LAYERS;	
	for (k=zStart;k<zEnd+1;k++){
		z = abs(zBound - k)%PML_LAYERS;
			
		p = PMLINDEX(i,j,k);
		dzStore = DZ(p);			
	
		m = MATERIALINDEX(i,j,k);
		eps = EPSR(m)*EPSNOT;
		DZ(p)     = cx[x].c1*DZ(p) + cx[x].c2*(HY(i,j,k) - HY(i-1,j,k) - HX(i,j,k) + HX(i,j-1,k) - JSZ(m));
		EZ(i,j,k) = cy[y].c3*EZ(i,j,k) + cy[y].c4*(cz[z].c5*DZ(p)-cz[z].c6*dzStore)/eps; 
	}}} // end for

}// end updatePMLEFieldUtilEz




void updatePMLEFields(int n){
	// figure out which pml layers are needed for the processor
	// All the conditions are independent of each other, for example, a west processor can also be a east processor!
	// so before each function call, we check every condition for that call
	int xStart, xEnd, yStart, yEnd, zStart, zEnd;

	// 6 faces

	// ux-WEST
	if (checkEdgeProcessor(WEST)){
		xStart = 2;
		xEnd   = PML_LAYERS;

		setStandardYCoordinates(&yStart, &yEnd);
		setStandardZCoordinates(&zStart, &zEnd);
		
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, 1, 1, cn, cOnes, cOnes, n);
		// for Boundary conditions at x = 1
		updatePMLEFieldUtilEx(1, 1, yStart, yEnd, zStart, zEnd, 1, 1, 1, cn, cOnes, cOnes, n);
	}// end Ux-West

	//Ux-EAST
	if (checkEdgeProcessor(EAST)){
		xStart = xLen-PML_LAYERS+1;
		xEnd   = xLen;

		setStandardYCoordinates(&yStart, &yEnd);
		setStandardZCoordinates(&zStart, &zEnd);
		
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, xLen, 1, 1, cn, cOnes, cOnes, n);
		
	}// end Ux-East

	// Uy-south
	if (checkEdgeProcessor(SOUTH)){
	
		yStart = 2;
		yEnd   = PML_LAYERS;
	
		setStandardXCoordinates(&xStart, &xEnd);
		setStandardZCoordinates(&zStart, &zEnd);
		
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, 1, 1, cOnes, cn, cOnes, n);
		// for boundary conditiosn at y=1
		updatePMLEFieldUtilEy(xStart, xEnd, 1, 1, zStart, zEnd, 1, 1, 1, cOnes, cn, cOnes, n);
		
	}// end Uy-South

	// Uy-NORTH
	if (checkEdgeProcessor(NORTH)){
		yStart = yLen-PML_LAYERS+1;
		yEnd   = yLen;

		setStandardXCoordinates(&xStart, &xEnd);
		setStandardZCoordinates(&zStart, &zEnd);
		
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, yLen, 1, cOnes, cn, cOnes, n);
		
	}// end Uy-NORTH

	// Uz-Bottom
	if (checkEdgeProcessor(BOTTOM)){
	
		zStart = 2;
		zEnd   = PML_LAYERS;	

		setStandardXCoordinates(&xStart, &xEnd);
		setStandardYCoordinates(&yStart, &yEnd);
		
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, 1, 1, cOnes, cOnes, cn, n);
		// for boundary conditiosn at z=1
		updatePMLEFieldUtilEz(xStart, xEnd, yStart, yEnd, 1, 1, 1, 1, 1, cOnes, cOnes, cn, n);
		
	}// end Uz-Bottom

	// Uz-TOP
	if (checkEdgeProcessor(TOP)){
		zStart = zLen-PML_LAYERS+1;
		zEnd   = zLen;

		setStandardXCoordinates(&xStart, &xEnd);
		setStandardYCoordinates(&yStart, &yEnd);

		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, 1, zLen, cOnes, cOnes, cn, n);
	}// End Uz-Top

	// Now the edges
//////////////////////////////////XY//////////////////
	//Uxy-west-south
        if (checkEdgeProcessor(WEST)&&checkEdgeProcessor(SOUTH)){
		xStart = 2;
		xEnd   = PML_LAYERS;
	
		yStart = 2;
		yEnd   = PML_LAYERS;

		setStandardZCoordinates(&zStart, &zEnd);

		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, 1, 1, cn, cn, cOnes, n);

	}// end U-SW	
	
	//Uxy-WEST-NORTH
        if (checkEdgeProcessor(WEST)&&checkEdgeProcessor(NORTH)){
		xStart = 2;
		xEnd   = PML_LAYERS;
	
		yStart = yLen-PML_LAYERS+1;
		yEnd   = yLen;

		setStandardZCoordinates(&zStart, &zEnd);
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, yLen, 1, cn, cn, cOnes, n);
		// X = 1
		updatePMLEFieldUtilEx(1, 1, yStart, yEnd, zStart, zEnd, 1, yLen, 1, cn, cn, cOnes, n);
	}// end Uxy-WN
	
	//Uxy-EAST-SOUTH
        if (checkEdgeProcessor(EAST)&&checkEdgeProcessor(SOUTH)){
	
		xStart = xLen-PML_LAYERS+1;
		xEnd   = xLen;
	
		yStart = 2;
		yEnd   = PML_LAYERS;

		setStandardZCoordinates(&zStart, &zEnd);
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, xLen, 1, 1, cn, cn, cOnes, n);
		// y = 1
		updatePMLEFieldUtilEy(xStart, xEnd, 1, 1, zStart, zEnd, xLen, 1, 1, cn, cn, cOnes, n);

	}// end Uxy-ES
		
	
	//Uxy-East-NORTH		
        if (checkEdgeProcessor(EAST)&&checkEdgeProcessor(NORTH)){
		xStart = xLen-PML_LAYERS+1;
		xEnd   = xLen;
	
		yStart = yLen-PML_LAYERS+1;
		yEnd   = yLen;

		setStandardZCoordinates(&zStart, &zEnd);
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, xLen, yLen, 1, cn, cn, cOnes, n);
	} // end U-North-East
////////////////////////////////YZ//////////////////
	
	//Uyz-South-BOTTTOM
        if (checkEdgeProcessor(SOUTH)&&checkEdgeProcessor(BOTTOM)){
		yStart = 2;
		yEnd   = PML_LAYERS;
	
		zStart = 2;
		zEnd   = PML_LAYERS;

		setStandardXCoordinates(&xStart, &xEnd);
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, 1, 1, cOnes, cn, cn, n);

	}// end Uyz-SOUTH-BOTTOM	

	//Uyz-SOUTH-TOP
        if (checkEdgeProcessor(SOUTH)&&checkEdgeProcessor(TOP)){
		yStart = 2;
		yEnd   = PML_LAYERS;
	
		zStart = zLen-PML_LAYERS+1;
		zEnd   = zLen;

		setStandardXCoordinates(&xStart, &xEnd);

		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, 1, zLen, cOnes, cn, cn, n);
		// y = 1
		updatePMLEFieldUtilEy(xStart, xEnd, 1, 1, zStart, zEnd, 1, 1, zLen, cOnes, cn, cn, n);

	}// end Uyz-SOUTH-TOP
	
	//Uyz-NORTH-BOTTTOM
        if (checkEdgeProcessor(NORTH)&&checkEdgeProcessor(BOTTOM)){
		yStart = yLen-PML_LAYERS+1;
		yEnd   = yLen;
	
		zStart = 2;
		zEnd   = PML_LAYERS;

		setStandardXCoordinates(&xStart, &xEnd);
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, yLen, 1, cOnes, cn, cn, n);
		// z=1
		updatePMLEFieldUtilEz(xStart, xEnd, yStart, yEnd, 1, 1, 1, yLen, 1, cOnes, cn, cn, n);
	}// end U-NORTH-Botttom

	//Uyz-NORTH-TOP
        if (checkEdgeProcessor(NORTH)&&checkEdgeProcessor(TOP)){
		yStart = yLen-PML_LAYERS+1;
		yEnd   = yLen;
	
		zStart = zLen-PML_LAYERS+1;
		zEnd   = zLen;

		setStandardXCoordinates(&xStart, &xEnd);
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, yLen, zLen, cOnes, cn, cn, n);
	}// end Uyz-NORTH-TOP

//////////////////ZX////////////////

	//Uzx-BOTTOM-WEST
        if (checkEdgeProcessor(WEST)&&checkEdgeProcessor(BOTTOM)){
		zStart = 2;
		zEnd   = PML_LAYERS;

		xStart = 2;
		xEnd   = PML_LAYERS;

		setStandardYCoordinates(&yStart, &yEnd);
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, 1, 1, cn, cOnes, cn, n);
	}// end Uzx-BOTTOM-WEST

	//Uzx-BOTTOM-EAST
        if (checkEdgeProcessor(BOTTOM)&&checkEdgeProcessor(EAST)){

		zStart = 2;
		zEnd   = PML_LAYERS;
		
		xStart = xLen-PML_LAYERS+1;
		xEnd   = xLen;
	
		setStandardYCoordinates(&yStart, &yEnd);
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, xLen, 1, 1, cn, cOnes, cn, n);
		//z=1
		updatePMLEFieldUtilEz(xStart, xEnd, yStart, yEnd, 1, 1, xLen, 1, 1, cn, cOnes, cn, n);
	}// end Uzx-BOTTOM-EAST

	//Uzx-TOP-WEST
        if (checkEdgeProcessor(TOP)&&checkEdgeProcessor(WEST)){
		zStart = zLen-PML_LAYERS+1;
		zEnd   = zLen;

		xStart = 2;
		xEnd   = PML_LAYERS;
	
		setStandardYCoordinates(&yStart, &yEnd);
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, 1, zLen, cn, cOnes, cn, n);
		//x=1
		updatePMLEFieldUtilEx(1, 1, yStart, yEnd, zStart, zEnd, 1, 1, zLen, cn, cOnes, cn, n);
	}// end Uzx-TOP-WEST

	//Uzx-TOP-EAST
        if (checkEdgeProcessor(TOP)&&checkEdgeProcessor(EAST)){
		zStart = zLen-PML_LAYERS+1;
		zEnd   = zLen;

		xStart = xLen-PML_LAYERS+1;
		xEnd   = xLen;
	
		setStandardYCoordinates(&yStart, &yEnd);
		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, xLen, 1, zLen, cn, cOnes, cn, n);
	}// end Uzx-TOP-EAST
	
	// 8 corners
//////////////////xyz///////////////
	//Uxyz-WEST-SOUTH-BOTTOM
        if (checkEdgeProcessor(WEST)&&checkEdgeProcessor(SOUTH)&&checkEdgeProcessor(BOTTOM)){
		xStart = 2;
		xEnd   = PML_LAYERS;

		yStart = 2;
		yEnd   = PML_LAYERS;
			
		zStart = 2;
		zEnd   = PML_LAYERS;

		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, 1, 1, cn, cn, cn, n);
	}// end Uxyz-WEST-SOUTH-BOTTOM
	
	//Uxyz-WEST-SOUTH-TOP
        if (checkEdgeProcessor(WEST)&&checkEdgeProcessor(SOUTH)&&checkEdgeProcessor(TOP)){
		xStart = 2;
		xEnd   = PML_LAYERS;

		yStart = 2;
		yEnd   = PML_LAYERS;
			
		zStart = zLen-PML_LAYERS+1;
		zEnd   = zLen;

		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, 1, zLen, cn, cn, cn, n);
	}// end Uxyz-WEST-SOUTH-TOP

	//Uxyz-WEST-NORTH-BOTTOM
        if (checkEdgeProcessor(WEST)&&checkEdgeProcessor(NORTH)&&checkEdgeProcessor(BOTTOM)){
		xStart = 2;
		xEnd   = PML_LAYERS;

		yStart = yLen-PML_LAYERS+1;
		yEnd   = yLen;
			
		zStart = 2;
		zEnd   = PML_LAYERS;

		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, yLen, 1, cn, cn, cn, n);
	}// end Uxyz-WEST-NORTH-BOTTOM

	// Uxyz-WEST-NORTH-TOP
        if (checkEdgeProcessor(WEST)&&checkEdgeProcessor(NORTH)&&checkEdgeProcessor(TOP)){
		xStart = 2;
		xEnd   = PML_LAYERS;

		yStart = yLen-PML_LAYERS+1;
		yEnd   = yLen;
			
		zStart = zLen-PML_LAYERS+1;
		zEnd   = zLen;

		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, 1, yLen, zLen, cn, cn, cn, n);
		// x=1;
		updatePMLEFieldUtilEx(1, 1, yStart, yEnd, zStart, zEnd, 1, yLen, zLen, cn, cn, cn, n);
	}// end Uxyz-WEST-NORTH-TOP

	//Uxyz-EAST-SOUTH-BOTTOM
        if (checkEdgeProcessor(EAST)&&checkEdgeProcessor(SOUTH)&&checkEdgeProcessor(BOTTOM)){
		xStart = xLen-PML_LAYERS+1;
		xEnd   = xLen;

		yStart = 2;
		yEnd   = PML_LAYERS;
			
		zStart = 2;
		zEnd   = PML_LAYERS;

		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, xLen, 1, 1, cn, cn, cn, n);
	}// end Uxyz-EAST-SOUTH-BOTTOM

	//Uxyz-EAST-SOUTH-TOP
        if (checkEdgeProcessor(EAST)&&checkEdgeProcessor(SOUTH)&&checkEdgeProcessor(TOP)){
		xStart = xLen-PML_LAYERS+1;
		xEnd   = xLen;

		yStart = 2;
		yEnd   = PML_LAYERS;
			
		zStart = zLen-PML_LAYERS+1;
		zEnd   = zLen;

		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, xLen, 1, zLen, cn, cn, cn, n);
		// y =1
		updatePMLEFieldUtilEy(xStart, xEnd, 1, 1, zStart, zEnd, xLen, 1, zLen, cn, cn, cn, n);
	}// end Uxyz-EAST-SOUTH-TOP

	
	// Uxyz-EAST-NORTH-BOTTOM
        if (checkEdgeProcessor(EAST)&&checkEdgeProcessor(NORTH)&&checkEdgeProcessor(BOTTOM)){
		xStart = xLen-PML_LAYERS+1;
		xEnd   = xLen;

		yStart = yLen-PML_LAYERS+1;
		yEnd   = yLen;
			
		zStart = 2;
		zEnd   = PML_LAYERS;

		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, xLen, yLen, 1, cn, cn, cn, n);
		// z=1;
		updatePMLEFieldUtilEz(xStart, xEnd, yStart, yEnd, 1, 1, xLen, yLen, 1, cn, cn, cn, n);
	}// end Uxyz-EAST-NORTH-BOTTOM

	// Uxyz-EAST-NORTH-TOP
        if (checkEdgeProcessor(EAST)&&checkEdgeProcessor(NORTH)&&checkEdgeProcessor(TOP)){
		xStart = xLen-PML_LAYERS+1;
		xEnd   = xLen;

		yStart = yLen-PML_LAYERS+1;
		yEnd   = yLen;
			
		zStart = zLen-PML_LAYERS+1;
		zEnd   = zLen;

		updatePMLEFieldUtilAll(xStart, xEnd, yStart, yEnd, zStart, zEnd, xLen, yLen, zLen, cn, cn, cn, n);
	}// end U-EAST-NORTH-TOP
	//printf("processor # %d end PMLE fields ",sid);	
}// end updatePMLEFields
		
void setStandardXCoordinates(int *xStart, int *xEnd){
		if (checkEdgeProcessor(WEST))
			*xStart = PML_LAYERS+1;
		else
			*xStart = 1;

		if (checkEdgeProcessor(EAST))
			*xEnd = xLen-PML_LAYERS;
		else
			*xEnd = xLen;

}// end  setStandardXCoordinates

void setStandardYCoordinates(int *yStart, int *yEnd){
		if (checkEdgeProcessor(SOUTH))
			*yStart = PML_LAYERS+1;
		else
			*yStart = 1;
		if (checkEdgeProcessor(NORTH))
			*yEnd = yLen-PML_LAYERS;
		else
			*yEnd = yLen;

}// end  setStandardYCoordinates


void setStandardZCoordinates(int *zStart, int *zEnd){
		if (checkEdgeProcessor(BOTTOM))
			*zStart = PML_LAYERS+1;
		else
			*zStart = 1;
		if (checkEdgeProcessor(TOP))
			*zEnd = zLen-PML_LAYERS;
		else
			*zEnd = zLen;

}// end  setStandardZCoordinates

