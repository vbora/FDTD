#include"common.h"
#include"defs.h"
#include"externGlobal.h"
#include<math.h>
#include<limits.h>
#include"util.h"

void setNeighboursAndCommunicationConsts(){
nebW=nebE=nebS=nebN=nebB=nebT=-1;

if (! checkEdgeProcessor(WEST))
        nebW = (px-1)*vproc[1]*vproc[2] + py*vproc[2] + pz;

if (! checkEdgeProcessor(EAST))
        nebE = (px+1)*vproc[1]*vproc[2] + py*vproc[2] + pz;

if (! checkEdgeProcessor(SOUTH))
        nebS = px*vproc[1]*vproc[2] + (py-1)*vproc[2] + pz;

if (! checkEdgeProcessor(NORTH))
        nebN = px*vproc[1]*vproc[2] + (py+1)*vproc[2] + pz;

if (! checkEdgeProcessor(BOTTOM))
        nebB = px*vproc[1]*vproc[2] + py*vproc[2] + pz-1;

if (! checkEdgeProcessor(TOP))
        nebT = px*vproc[1]*vproc[2] + py*vproc[2] + pz+1;

copyCountEW = yLen * zLen * 2;  
copyCountNS = xLen * zLen * 2; 
copyCountTB = xLen * yLen * 2;


}// neighbours

void errorFunction(char a[100]){
	// TODO print more things like sid etc
	printf("\nError occured --> %s processor # %d\n",a, sid);
	exit (-1);
} // End error function

int** alloc2DArray(size_t size_x, size_t size_y, size_t elemSize){
	int **a;
	int i;
	a = calloc(size_x, sizeof(int*));
	if (a == NULL){
		printf("Before 1st error function sid %d", sid);
		errorFunction("Malloc failed");		
	}
	for (i=0; i < size_x; i++){
		a[i] =  calloc(size_y , elemSize);
		if (a[i] == NULL){
			printf("Before 2nd error function sid %d", sid);
			errorFunction("Malloc failed");		
			}
	}
	return a; 

}// end alloc2DDoubleArray 

int*** alloc3DArray(size_t sizeX, size_t sizeY, size_t sizeZ, size_t elemSize){
	int ***a;
	int i,j;
	a = calloc(sizeX, sizeof(int**));
	if (a == NULL){
		printf("Before 1st (3d) error function sid %d size %d \n", sid, sizeX);
		errorFunction("Malloc failed initially");		
		}
	for (i=0; i < sizeX; i++){
		a[i] =  calloc(sizeY , sizeof(int*));
		if (a[i] == NULL){
			printf("Before 2nd (3d) error function sid %d", sid);	
			errorFunction("Malloc failed second step");	
		}
		for (j=0; j < sizeY; j++){
			a[i][j] =  calloc(sizeZ , elemSize);
			if (a[i][j] == NULL){
				printf("Before 3rd (3d) error function sid %d", sid);	
				errorFunction("Malloc failed final allocation");	
			}
		
		}	
	}
	return a; 

}// end alloc3DArray 

void maxInt(int x, int y, int z, int *largest, int *secondLargest, int *smallest ){
	if (x > y){
		if (x>z){
			*largest = x;
			if (y>z){
				*secondLargest = y;
				*smallest = z;
			}
			else{
				*secondLargest = z;
				*smallest = y;
			}
		}
	}
	else{
		if (y>z){
			*largest = y;
			if (x>z){
				*secondLargest = x;
				*smallest = z;
			}
			else{
				*secondLargest = z;
				*smallest = x;
			}
		}
		else{ 
			*largest = z;
			if (x>y){
				*secondLargest = x;
				*smallest = y;
			}
			else{
				*secondLargest = y;
				*smallest = x;
			}
		}
	}
}// end function

int isLocal(int gx, int gy, int gz, int *ix, int *iy, int *iz){
	// The global coordinates also start from (0,0,0), however the local coordinates start from 1,1 (Since there is an extra row for the communication
	// So we add the 1 while calculating local
	if (( gx / xLen == px ) && (gy / yLen == py) && (gz / zLen ==  pz) ){
	        *ix = gx % xLen + 1; 
		*iy = gy % yLen + 1;
		*iz = gz % zLen + 1;	
		return 1;
	}
	else{
		*ix = -1; 
		*iy = -1;
		*iz = -1;
		return 0;
	}
} //end isLocal


void globalEquivalent(int ix, int iy, int iz, int *gx, int *gy, int *gz){
	// Given local coordinates, find the global coordinates
	// The global coordinates also start from (0,0,0), however the local coordinates start from 1,1 (Since there is an extra row for the communication
	*gx = px*xLen + ix-1;
	*gy = py*yLen + iy-1;
	*gz = pz*zLen + iz-1;

} //end globalEquivalent


void makeSphere(int m, int center[3], int radius){

	int gx,gy,gz,ix,iy,iz;
	int xStart, xEnd, yStart, yEnd, zStart, zEnd;
	int rSquare;
	int sphNumber;
	// drawing sphere 1 of material m
/*	xStart = 0; 
	xEnd = XDIM-1; 
	yStart = 0;
	yEnd = YDIM -1;  */

	xStart = center[X] - radius - 1 ;// since index from 0
	xEnd   = center[X] + radius + 1; // so that I can use <
	yStart = center[Y] - radius - 1 ;// since index from 0
	yEnd   = center[Y] + radius + 1; // so that I can use <
	zStart = center[Z] - radius - 1 ;// since index from 0
	zEnd   = center[Z] + radius + 1; // so that I can use <
	
	//printf("xStart %d, xEnd %d  yStart %d yEnd %d zStart %d zEnd %d\n",xStart, xEnd, yStart, yEnd, zStart, zEnd);	
	//printf( "center(%d,%d,%d)\n", center[X], center[Y], center[Z]);
	if (center[Y] == 10)	
		sphNumber = 1;
	else
		sphNumber = 2;
	rSquare = radius*radius;	
	for (gx=xStart; gx <= xEnd; gx++)
	for (gy=yStart; gy <= yEnd; gy++)
	for (gz=zStart; gz <= zEnd; gz++){
		if (rSquare >= ((gx-center[X])*(gx-center[X]) + (gy-center[Y])*(gy-center[Y])+(gz-center[Z])*(gz-center[Z]))){
			if (isLocal(gx,gy,gz,&ix,&iy,&iz)){

				MATERIALINDEX(ix,iy,iz)=m;
				/*
				if (NPOLES(m) !=0)
					DISPFIELDINDEX(ix,iy,iz)=dispFieldIndex++;
				*/
			}
		}
	} // end for
}// end make Sphere 


// reads global coordinates and vproc
// To check if the processor is at the "side" edge of the simulation space
int checkEdgeProcessor(int side){
/*

	North
	 __	
West  	|  | East
  	|__|
  (0,0)	South

an TOP BOTTOM(z=zero) on the z	
 
*/

	int normVProc[3];
	normVProc[0] = vproc[0] - 1;
	normVProc[1] = vproc[1] - 1;
	normVProc[2] = vproc[2] - 1;	
	switch(side){
		case NORTH: 
			if (py == normVProc[1])
				return 1;
			else 
				return 0;
		case SOUTH: 
			if (py == 0)
				return 1;
			else 
				return 0;
		case WEST: 
			if (px == 0)
				return 1;
			else 
				return 0;
		     	
		case EAST: 
			if (px == normVProc[0])
				return 1;
			else 
				return 0;
	
		case BOTTOM:
			if (pz == 0)
				return 1;
			else 
				return 0;
		case TOP:
			if (pz == normVProc[2])
				return 1;
			else 
				return 0;
				
	}// end switch-case
	return 0;
} 



void splitProb(int nProc, int xDim, int yDim, int zDim, int finProc[3]) {
	int i,j,num;
	int subFactors[MAX_FACT];	
	int subNumberofFactors;
	int n1,n2,n3;
	int subNumber;
	double sideX,sideY,sideZ;
	int numberOfFactors,factors[MAX_FACT];
	// TODO make this ULLONG_MAX
	unsigned long long comm = ULONG_MAX;
	unsigned long long commTest;

	numberOfFactors = factorize(nProc, factors);

	for (i = 0; i < numberOfFactors; i++){

		n1 = factors[i];	
		subNumber = nProc / n1;		
		subNumberofFactors = factorize(subNumber, subFactors);

		for (j = 0; j < subNumberofFactors; j++){
			n2 = subFactors[j];
			n3 = subNumber / n2;

			// checking if grid points > number of processors
			
			if ((xDim < n1) || (yDim < n2) || (zDim < n3)){
				//printf("\n grid points smaller than processors %d %d %d", n1,n2,n3);	
				continue;
			}			
			sideX =	(double) xDim / n1;
			sideY =	(double) yDim / n2;
			sideZ =	(double) zDim / n3;

			// checking for integer values
			if ((sideX > (int)sideX) || (sideY > (int)sideY) || (sideZ > (int)sideZ)){
				//printf("\n not integer grid points %d %d %d", n1,n2,n3);	
				continue;
			}
			else{ 
				
				commTest = sideX*sideX + sideY*sideY + sideZ*sideZ;				
				//printf("\n comm test %lu", commTest);
				//printf("\n %d %d %d", n1,n2,n3);	
				if (commTest < comm){
					comm = commTest;
					finProc[0] = n1;
					finProc[1] = n2;
					finProc[2] = n3;
				}

			}
		}
			

	}// end for
	if (sid == 0)
		printf("problem split into %d %d %d", finProc[0], finProc[1], finProc[2]);	
}// end splitProb

int factorize(int inNum, int factors[MAX_FACT]){
int i, j,num, numberOfFactors ;
double sqrtNum;
num = inNum;	


numberOfFactors = 0;
i = 1;
	while (i <= num)
	{	
		if (num % i == 0){
			factors[numberOfFactors]=i;
			//num = num / i;
			numberOfFactors++;
		}
		i++;
			
	}// end while
	
	return numberOfFactors;
}

