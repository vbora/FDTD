#include"common.h"
#include"defs.h"
#include"globals.h"
#include"main.h"
#include"fileUtil.h"
#include"util.h"

int main(int argc, char **argv){
//	// E and H are already intialized to zero, 
	// Hence only initialize Ca and Da
        time_t startTime ,endTime;
        time_t tempStartTime , tempEndTime;
        double runTime;
	int n, totalTimeSteps = TOTALTIMESTEPS;
	char filePrefix[100];
	int m;
        // for monitor points
        int i,j,k;
        int gx,gy,gz;


	MPI_Status sendStatus1, sendStatus2, recvStatus1, recvStatus2;
	MPI_Request recvRequest1, recvRequest2; 
 	
	// There are 2 sets of coordinates, local (internal) coordinates, and global coordinates
	// Namesly gx, gy, gz, and local as ix, iy, iz 
	
	// The processors are defined by px, py, pz
	// The relationships are
	// The dimensions of the big simulation box (Global martix XDIM, YDIM, ZDIM
	// The dimension of the grid in the processors are iXdim, iYdim, iZdim
 
	// gx = px * iXdim
	// gy = py * iYdim
	// gz = pz * iZdim

	gx = XDIM/2-1;
        gy = YDIM/2-1;
        gz = ZDIM/2-1;


        time(&startTime);
	init(argc, argv);

	writeMatFile("MAT");
	for (n = 0; n < totalTimeSteps; n++){
		//printf("\nTime step %d processor # %d", n,sid);
			
		updateSource(n);
		//printf("\nTime step %d after update source processor # %d", n,sid);
		time(&tempStartTime);

		processStep(n);
		time(&tempEndTime);
        	runTime = difftime(tempEndTime,tempStartTime);
        	//printf("\n processor # %d STEP NUMBER %d step time %e",sid, n, runTime);
		//printStep(n);
		 if (isLocal(gx, gy, gz, &i, &j, &k))
                        printf("\n TimeStep %d, value Ex %e Ey %e", n, EX(i,j,k),EY(i,j,k));

	}// end for
        time(&endTime);
        runTime = difftime(endTime,startTime);
        printf("\n processor # %d TOTAL TIME %e",sid, runTime);
	//free your memory, do you have to free each pointer?
	MPI_Finalize();
	// have a function to free up everything
		
}// end main
