#include"common.h"
#include"defs.h"
#include"externGlobal.h"
#include"util.h"
#include"updateFields.h" 
#include"updatePMLFields.h" 
#include"dispersiveMat.h"
void processStep(int n){
	time_t startTime ,endTime;
	double runTime;
	//printf("\nPROCESSOR %d TIMESTEP %d inside process step", sid, n);	
	time(&startTime);
	updateEFields(n);
	time(&endTime);
	runTime = difftime(endTime,startTime);

	//printf("\nPROCESSOR %d TIMESTEP %d update E fields time %e", sid, n, runTime);

	time(&startTime);
	if (PML_LAYERS > 0)
		updatePMLEFields(n);
	time(&endTime);
	runTime = difftime(endTime,startTime);

	//printf("\nPROCESSOR %d TIMESTEP %d PML E fields time %e", sid, n,  runTime);
	time(&startTime);
	communicateEfield(n);
	time(&endTime);
	runTime = difftime(endTime,startTime);

	//printf("\nPROCESSOR %d TIMESTEP %d communicate E fields time %e", sid,n, runTime);

	
	updateHFields(n);
	//printf("\nsid %d  TIMESTEP %d after updateHFields nebE %d ",sid,n, nebE);	
	if (PML_LAYERS > 0)
		updatePMLHFields(n);
	//printf("\nsid %d TIMESTEP %d after updatePMLHFields nebE %d ",sid,n, nebE);	
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
	int xEnd, yEnd, zEnd;
	// In global coordinates 
	// We donot write Ex[0][j] and Ey[i][0], so as to keep them zero 
	// Hence takes care of all sides and all corners of Ex and Ey
	
	xStart = yStart = zStart = 1 ;
	xEnd = xLen; 
	yEnd = yLen; 
	zEnd = zLen;

	if(PML_LAYERS==0){
		if(checkEdgeProcessor(WEST))
			xStart = 2;
		if (checkEdgeProcessor(SOUTH)) 
			yStart = 2;
		if(checkEdgeProcessor(BOTTOM))
			zStart = 2;
	
		if(checkEdgeProcessor(EAST))
			xEnd = xLen;
		if (checkEdgeProcessor(NORTH)) 
			yEnd = yLen;
		if(checkEdgeProcessor(TOP))
			zEnd = zLen;
	}
	else{
		if(checkEdgeProcessor(WEST))
			xStart = 1+PML_LAYERS;
		if (checkEdgeProcessor(SOUTH)) 
			yStart = 1+PML_LAYERS;
		if(checkEdgeProcessor(BOTTOM))
			zStart = 1+PML_LAYERS;
	
		if(checkEdgeProcessor(EAST))
			xEnd = xLen-PML_LAYERS; 
		if(checkEdgeProcessor(NORTH))
			yEnd = yLen-PML_LAYERS;
		if(checkEdgeProcessor(TOP))
			zEnd = zLen-PML_LAYERS;
	} // end if

	//printf("\nPROCESSOR # %d xStart %d, xEnd %d,  yStart %d, yEnd %d,  zStart %d, zEnd %d",sid, xStart ,xEnd ,yStart, yEnd,  zStart , zEnd);
	// The total Ex/Ey/Hz are till Ex[xLen+1][yLen+1] but the last rows are not this processors, hence read only 	
	for( i=xStart;i<xEnd+1; i++)
	for( j=yStart;j<yEnd+1; j++)
	for( k=zStart;k<zEnd+1; k++){
		m = MATERIALINDEX(i,j,k);
		if (NPOLES(m)==0){
			EX(i,j,k) = CA(m) * EX(i,j,k) + CB(m)*(HZ(i,j,k) - HZ(i,j-1,k) + HY(i,j,k-1) - HY(i,j,k) - JSX(m) );		  
			EY(i,j,k) = CA(m) * EY(i,j,k) + CB(m)*(HX(i,j,k) - HX(i,j,k-1) + HZ(i-1,j,k) - HZ(i,j,k) - JSY(m) );		 
			EZ(i,j,k) = CA(m) * EZ(i,j,k) + CB(m)*(HY(i,j,k) - HY(i-1,j,k) + HX(i,j-1,k) - HX(i,j,k) - JSZ(m) ); 
			if (EY(i,j,k) > maxE)
				maxE = EY(i,j,k);
			//if (EX(i,j,k) != 0)
			//	printf("TIME STEP %d value %e, (%d,%d,%d) ", n, EX(i,j,k), i,j,k);		maxE = EX(i,j,k);
		}
		else 
			updateDispField(i,j,k,m,n);
	}// end for



	// PLANES
	//The end planes, TOP, NORTH and EAST are already done

	// The tangential components are zero, and we need to calculate only the normal components
	if (PML_LAYERS==0){
		// BOTTOM
		if (checkEdgeProcessor(BOTTOM)){ 
			k = 1 ; //--> Boundary at BOTTOM (z = 1) 
			for( i=xStart; i<xEnd+1; i++)
			for( j=yStart; j<yEnd+1; j++){
				m = MATERIALINDEX(i,j,k);
				EZ(i,j,k) = CA(m) * EZ(i,j,k) + CB(m)*(HY(i,j,k) - HY(i-1,j,k) + HX(i,j-1,k) - HX(i,j,k) - JSZ(m) ); 
			}
		}
		// WEST
		if (checkEdgeProcessor(WEST)){ 
			i = 1 ; //--> Boundary at WEST ( x = 1 )
			for( j=yStart; j<yEnd+1; j++)
			for( k=zStart; k<zEnd+1; k++){
				m = MATERIALINDEX(i,j,k);
				EX(i,j,k) = CA(m) * EX(i,j,k) + CB(m)*(HZ(i,j,k) - HZ(i,j-1,k) + HY(i,j,k-1) - HY(i,j,k) - JSX(m) );		  
			}
		}
		// SOUTH
		if (checkEdgeProcessor(SOUTH)){ 
			j = 1 ; //--> Boundary at SOUTH ( y = 1 )
			for( i=xStart; i<xEnd+1; i++)
			for( k=zStart; k<zEnd+1; k++){
				m = MATERIALINDEX(i,j,k);
				EY(i,j,k) = CA(m) * EY(i,j,k) + CB(m)*(HX(i,j,k) - HX(i,j,k-1) + HZ(i-1,j,k) - HZ(i,j,k) - JSY(m) );		 
				if (EY(i,j,k) > maxE)
					maxE = EY(i,j,k);
			}
		}
	
		//printf("\tprocessor # %d, Ey field %e\n", sid,maxE);

		// Corners
	} // End PML == 0

	
	//printf("\t processor # %d timestep %d max E fields %e\n", sid, n, maxE);
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
	int xEnd, yEnd, zEnd;

	maxH = minH = 0;
	// because we use < not <= and the xLen+1 element is from a diff processor
	xStart = yStart = zStart =  1;	
	
	xEnd = xLen; 
	yEnd = yLen; 
	zEnd = zLen;

	if(PML_LAYERS==0){
		if(checkEdgeProcessor(WEST))
			xStart = 2;
		if (checkEdgeProcessor(SOUTH)) 
			yStart = 2;
		if(checkEdgeProcessor(BOTTOM))
			zStart = 2;
	
		if (checkEdgeProcessor(NORTH)) 
			yEnd = yLen;
		if(checkEdgeProcessor(EAST))
			xEnd = xLen;
		if(checkEdgeProcessor(TOP))
			zEnd = zLen;
	}
	else{
		if(checkEdgeProcessor(WEST))
			xStart = 1+PML_LAYERS;
		if (checkEdgeProcessor(SOUTH)) 
			yStart = 1+PML_LAYERS;
		if(checkEdgeProcessor(BOTTOM))
			zStart = 1+PML_LAYERS;
	
		if(checkEdgeProcessor(EAST))
			xEnd = xLen-PML_LAYERS; 
		if (checkEdgeProcessor(NORTH))
			yEnd = yLen-PML_LAYERS;
		if(checkEdgeProcessor(TOP))
			zEnd = zLen-PML_LAYERS;
	} // end if

	// Processes all points for non boundary processors
	for( i=xStart; i<xEnd+1; i++)
	for( j=yStart; j<yEnd+1; j++)
	for( k=zStart; k<zEnd+1; k++){
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

	if (PML_LAYERS==0){
		if (checkEdgeProcessor(BOTTOM)){ 
			//maxH = 0;
			k = 1;
			for (i=xStart ; i < xEnd+1; i++)
			for (j=yStart ; j < yEnd+1; j++){
				m = MATERIALINDEX(i,j,k);
				HX(i,j,k) = DA(m) * HX(i,j,k) + DB(m) * (EY(i,j,k+1) - EY(i,j,k) + EZ(i,j,k) - EZ(i,j+1,k) - MSX(m) );	
				HY(i,j,k) = DA(m) * HY(i,j,k) + DB(m) * (EZ(i+1,j,k) - EZ(i,j,k) + EX(i,j,k) - EX(i,j,k+1) - MSY(m) );
			}
		}

		if(checkEdgeProcessor(SOUTH)){
			//maxH = 0;
			j = 1;
			for (i=xStart; i < xEnd+1; i++)
			for (k=zStart; k < zEnd+1; k++){
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
			for (j=yStart; j < yEnd+1; j++)
			for (k=zStart; k < zEnd+1; k++){
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
			for (i=xStart; i < xEnd+1; i++){
				m = MATERIALINDEX(i,j,k);
				HX(i,j,k) = DA(m) * HX(i,j,k) + DB(m) * (EY(i,j,k+1) - EY(i,j,k) + EZ(i,j,k) - EZ(i,j+1,k) - MSX(m) );	
			}
		}	
		if(checkEdgeProcessor(WEST)&& checkEdgeProcessor(BOTTOM)){
			i=1;
			k=1;
			for (j=yStart; j < yEnd+1; j++){
				m = MATERIALINDEX(i,j,k);
				HY(i,j,k) = DA(m) * HY(i,j,k) + DB(m) * (EZ(i+1,j,k) - EZ(i,j,k) + EX(i,j,k) - EX(i,j,k+1) - MSY(m) );
			}
		}
	
		if(checkEdgeProcessor(WEST)&& checkEdgeProcessor(SOUTH)){
			i=1;
			j=1;
			for (k=zStart; k < zEnd+1; k++){
				m = MATERIALINDEX(i,j,k);
				HZ(i,j,k) = DA(m) * HZ(i,j,k) + DB(m) * (EX(i,j+1,k) - EX(i,j,k) + EY(i,j,k) - EY(i+1,j,k) - MSZ(m) );		  
				if (HZ(i,j,k) > maxH)
					maxH = HZ(i,j,k);
			}
		}		
		
	
		//if(checkEdgeProcessor(WEST)&& checkEdgeProcessor(SOUTH)&& checkEdgeProcessor(BOTTOM))
		// Everything is zero	
	}// End PML == 0
	//printf("\t processor # %d timestep %d max hFields %e\n", sid, n, maxH);
}// end updateHFields
