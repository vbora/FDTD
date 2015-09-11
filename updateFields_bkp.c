#include"common.h"
#include"defs.h"
#include"externGlobal.h"
#include"util.h"
#include"updateFields.h" 

void processStep(int n){
	updateEFields(n);
	communicateEfield(n);
	updateHFields(n);
	communicateHfield(n);

/* edge3  side3 edge2	
	|-------| 
 side4	|	| side 2
	|	|
	|-------|
  edge4	  side1  edge1

Our unit cell
	 - - - 
      ^	|     |
    Ey|	|  OHz|	
	|     |	
	 - - - 
	  ->Ex	

 Hence we just have to make the tangential Ex, Ey = 0 at the boundary (Since Hz is inside)
*/

}// processStep

void updateEFields(int n){
	int i,j,k,m;
	double maxE=0;
	int xStart, yStart, zStart;
	// In global coordinates 
	// We donot write Ex[0][j] and Ey[i][0], so as to keep them zero 
	// Hence takes care of all sides and all corners of Ex and Ey
	xStart = yStart = zStart = 1;

	if (checkEdgeProcessor(SOUTH)) 
		yStart = 2;
	if(checkEdgeProcessor(WEST))
		xStart = 2;
	if(checkEdgeProcessor(BOTTOM))
		zStart = 2;
		 	
	// The total Ex/Ey/Hz are till Ex[xLen+1][yLen+1] but the last rows are not this processors, hence read only 	
	for( i=xStart;i<xLen+1; i++)
	for( j=yStart;j<yLen+1; j++)
	for( k=zStart;k<zLen+1; k++){
		m = MATERIALINDEX(i,j,k);

		EX(i,j,k) = CA(m) * EX(i,j,k) + CB(m)*(HZ(i,j,k) - HZ(i,j-1,k) - HY(i,j,k) + HY(i,j,k-1) - JSX(m) );		  
		EY(i,j,k) = CA(m) * EY(i,j,k) + CB(m)*(HX(i,j,k) - HX(i,j,k-1) + HZ(i-1,j,k) - HZ(i,j,k) - JSY(m) );		 
		EZ(i,j,k) = CA(m) * EZ(i,j,k) + CB(m)*(HY(i,j,k) - HY(i-1,j,k) - HX(i,j,k) + HX(i,j-1,k) - JSZ(m) ); 
		if (EY(i,j,k) > maxE)
			maxE = EY(i,j,k);
		//if (EX(i,j,k) != 0)
		//	printf("TIME STEP %d value %e, (%d,%d,%d) ", n, EX(i,j,k), i,j,k);		maxE = EX(i,j,k);
	}// end for

	// PLANES
	//The end planes, TOP, NORTH and EAST are already done

	// The tangential components are zero, and we need to calculate only the normal components

	// BOTTOM
	if (checkEdgeProcessor(BOTTOM)){ 
		k = 1 ; //--> Boundary at BOTTOM (z = 1) 
		for( i=xStart; i<xLen+1; i++)
		for( j=yStart; j<yLen+1; j++){
			m = MATERIALINDEX(i,j,k);
			EZ(i,j,k) = CA(m) * EZ(i,j,k) + CB(m)*(HY(i,j,k) - HY(i-1,j,k) - HX(i,j,k) + HX(i,j-1,k) - JSZ(m) ); 
		}
	}
	// WEST
	if (checkEdgeProcessor(WEST)){ 
		i = 1 ; //--> Boundary at WEST ( x = 1 )
		for( j=yStart; j<yLen+1; j++)
		for( k=zStart; k<zLen+1; k++){
			m = MATERIALINDEX(i,j,k);
			EX(i,j,k) = CA(m) * EX(i,j,k) + CB(m)*(HZ(i,j,k) - HZ(i,j-1,k) - HY(i,j,k) + HY(i,j,k-1) - JSX(m) );		  
		}
	}
	// SOUTH
	if (checkEdgeProcessor(SOUTH)){ 
		j = 1 ; //--> Boundary at SOUTH ( y = 1 )
		for( i=xStart; i<xLen+1; i++)
		for( k=zStart; k<zLen+1; k++){
			m = MATERIALINDEX(i,j,k);
			EY(i,j,k) = CA(m) * EY(i,j,k) + CB(m)*(HX(i,j,k) - HX(i,j,k-1) + HZ(i-1,j,k) - HZ(i,j,k) - JSY(m) );		 
			if (EY(i,j,k) > maxE)
				maxE = EY(i,j,k);
		}
	}
	
	printf("\t processor # %d timestep %d max E fields %e\n", sid, n, maxE);
	//printf("\tprocessor # %d, Ey field %e\n", sid,maxE);

// Corners


	
}// end updateEFields

void updateHFields(int n){
/*        North  	
	|-------| 
 West	|	| East
	|	|
	|-------|
  	  south  
Our unit cell
*/
	int i,j,k,m; 
	double maxH, minH;
	int xStart, yStart, zStart; 

	maxH = minH = 0;
	// because we use < not <= and the xLen+1 element is from a diff processor
	xStart = 1;	
	yStart = 1;	
	zStart = 1;	
	
	if (checkEdgeProcessor(SOUTH)) 
		yStart  = 2 ;
	if(checkEdgeProcessor(WEST))
		xStart = 2 ;
	if(checkEdgeProcessor(BOTTOM))
		zStart = 2 ;
	// Processes all points for non boundary processors
	for( i=xStart; i<xLen+1; i++)
	for( j=yStart; j<yLen+1; j++)
	for( k=zStart; k<zLen+1; k++){
		m = MATERIALINDEX(i,j,k);
		HX(i,j,k) = DA(m) * HX(i,j,k) + DB(m) * (EY(i,j,k+1) - EY(i,j,k) + EZ(i,j,k) - EZ(i,j+1,k) - MSX(m) );	
		HY(i,j,k) = DA(m) * HY(i,j,k) + DB(m) * (EZ(i+1,j,k) - EZ(i,j,k) + EX(i,j,k) - EX(i,j,k+1) - MSY(m) );
		HZ(i,j,k) = DA(m) * HZ(i,j,k) + DB(m) * (EX(i,j+1,k) - EX(i,j,k) + EY(i,j,k) - EY(i+1,j,k) - MSZ(m) );
		//printf("HZ(i,j,k) %e,  DA(m) %e, DB(m) %e MSZ(m) %e \n",HZ(i,j,k), DA(m), DB(m), MSZ(m)); 
		if (HZ(i,j,k) > maxH)
			maxH = HZ(i,j,k);
		if (HZ(i,j,k) < minH)
			minH = HZ(i,j,k);
	}// end for
       	//printf("\t processor # %d timestep  %d minimum h  %e \n", sid, n, minH);
	//printf("\t processor # %d hz field in main %e \n", sid, maxH);

	if (checkEdgeProcessor(BOTTOM)){ 
		//maxH = 0;
		k = 1;
		for (i=xStart ; i < xLen+1; i++)
		for (j=yStart ; j < yLen+1; j++){
			m = MATERIALINDEX(i,j,k);
			HX(i,j,k) = DA(m) * HX(i,j,k) + DB(m) * (EY(i,j,k+1) - EY(i,j,k) + EZ(i,j,k) - EZ(i,j+1,k) - MSX(m) );	
			HY(i,j,k) = DA(m) * HY(i,j,k) + DB(m) * (EZ(i+1,j,k) - EZ(i,j,k) + EX(i,j,k) - EX(i,j,k+1) - MSY(m) );
		}
	}

	if(checkEdgeProcessor(SOUTH)){
		//maxH = 0;
		j = 1;
		for (i=xStart; i < xLen+1; i++)
		for (k=zStart; k < zLen+1; k++){
			m = MATERIALINDEX(i,j,k);
			HX(i,j,k) = DA(m) * HX(i,j,k) + DB(m) * (EY(i,j,k+1) - EY(i,j,k) + EZ(i,j,k) - EZ(i,j+1,k) - MSX(m) );	
			HZ(i,j,k) = DA(m) * HZ(i,j,k) + DB(m) * (EX(i,j+1,k) - EX(i,j,k) + EY(i,j,k) - EY(i+1,j,k) - MSZ(m) );		  
			if (HZ(i,j,k) > maxH)
				maxH = HZ(i,j,k);
		}
		//printf("\t processor # %d hz field in SOUTH %e \n", sid, maxH);
	}
	if(checkEdgeProcessor(WEST)){
		//maxH = 0;
		i = 1;
		for (j=yStart; j < yLen+1; j++)
		for (k=zStart; k < zLen+1; k++){
			m = MATERIALINDEX(i,j,k);
			HY(i,j,k) = DA(m) * HY(i,j,k) + DB(m) * (EZ(i+1,j,k) - EZ(i,j,k) + EX(i,j,k) - EX(i,j,k+1) - MSY(m) );
			HZ(i,j,k) = DA(m) * HZ(i,j,k) + DB(m) * (EX(i,j+1,k) - EX(i,j,k) + EY(i,j,k) - EY(i+1,j,k) - MSZ(m) );		  
			if (HZ(i,j,k) > maxH)
				maxH = HZ(i,j,k);
		}
		//printf("\t processor # %d hz field in WEST %e \n", sid, maxH);
	}
	// EDGES
	if(checkEdgeProcessor(SOUTH)&& checkEdgeProcessor(BOTTOM)){
		j=1;
		k=1;
		for (i=xStart; i < xLen+1; i++){
			m = MATERIALINDEX(i,j,k);
			HX(i,j,k) = DA(m) * HX(i,j,k) + DB(m) * (EY(i,j,k+1) - EY(i,j,k) + EZ(i,j,k) - EZ(i,j+1,k) - MSX(m) );	
		}
	}	
	if(checkEdgeProcessor(WEST)&& checkEdgeProcessor(BOTTOM)){
		i=1;
		k=1;
		for (j=yStart; j < yLen+1; j++){
			m = MATERIALINDEX(i,j,k);
			HY(i,j,k) = DA(m) * HY(i,j,k) + DB(m) * (EZ(i+1,j,k) - EZ(i,j,k) + EX(i,j,k) - EX(i,j,k+1) - MSY(m) );
		}
	}
	
	if(checkEdgeProcessor(WEST)&& checkEdgeProcessor(SOUTH)){
		i=1;
		j=1;
		for (k=zStart; k < zLen+1; k++){
			m = MATERIALINDEX(i,j,k);
			HZ(i,j,k) = DA(m) * HZ(i,j,k) + DB(m) * (EX(i,j+1,k) - EX(i,j,k) + EY(i,j,k) - EY(i+1,j,k) - MSZ(m) );		  
			if (HZ(i,j,k) > maxH)
				maxH = HZ(i,j,k);
		}
	}	
		
	
	//if(checkEdgeProcessor(WEST)&& checkEdgeProcessor(SOUTH)&& checkEdgeProcessor(BOTTOM))
	// Everything is zero	
	printf("\t processor # %d timestep %d max hFields %e\n", sid, n, maxH);
	
}// end updateHFields
