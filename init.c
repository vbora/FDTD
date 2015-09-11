#include"common.h"
#include"defs.h"
#include"externGlobal.h"
#include"util.h"
#include"init.h"
#include"dispersiveMat.h"

void memoryAllocation(size_t sizeX, size_t sizeY, size_t sizeZ){
	int largest, secondLargest, smallest;
	int bufSize;
	fieldMat = (fieldMaterialStruct***)alloc3DArray(sizeX,sizeY,sizeZ, sizeof(fieldMaterialStruct));  	
	printf("Memory allocated for field processor %d\n",sid);
	
	//allocate memory for the send and recv buffers
	maxInt(sizeX,sizeY,sizeZ, &largest, &secondLargest, &smallest);	
	
	bufSize = 6 * (largest) * (secondLargest) ;
	//printf ("\n sizeX %d, sizeY %d, sizeZ %d, largest %d, secondLargest %d , bufSize %d\n",sizeX, sizeY, sizeZ, largest, secondLargest, bufSize);

	sendBuffer1 = (double*)malloc(bufSize* sizeof(double));
	recvBuffer1 = (double*)malloc(bufSize* sizeof(double));

	sendBuffer2 = (double*)malloc(bufSize* sizeof(double));
	recvBuffer2 = (double*)malloc(bufSize* sizeof(double));
	
	
	//printf("In memory allocation, after all the buffers processor $ %d\n",sid);


	if (sendBuffer1 == NULL || recvBuffer1 == NULL || sendBuffer2 == NULL || recvBuffer2 == NULL)
		errorFunction("Malloc failed");
				
} // end memoryAllocation

void init(int argc, char ** argv){
	// set geometry constraints / read from file here
	// Do all the eps averaging, difficult since using the eps averaging thingy, so m increases a lot
        time_t startTime ,endTime;
        double runTime;	
	int nProc;		
	int i,j;
	int bufSize;
	int wrap_around[3]; // Flag to MPI_Cart_Create to tell it if our boundaries are periodic
	int reorder = 0; // If reordering is okay, Check with TODO Adam Mock
	
	
		
	
	wrap_around[0] = wrap_around[1] =wrap_around[2] = 0;

        MPI_Init(&argc,&argv); ///Initialize the MPI environment 
	MPI_Comm_size(MPI_COMM_WORLD, &nProc );
	splitProb(nProc,XDIM, YDIM, ZDIM, vproc);

	
	MPI_Cart_create(MPI_COMM_WORLD,3,vproc,wrap_around, reorder, &grid_comm);
	MPI_Comm_rank(grid_comm, &sid);	
	MPI_Cart_coords(grid_comm, sid, 3, coordinates);
	
	px = coordinates[0];
	py = coordinates[1];
	pz = coordinates[2];

	printf(" my sid %d , my coordinates %d , %d %d:\n",sid, px,py,pz);
	
	// will fail for 2D	
	deltaT = S * DELTA / c;	// S --> Courant factor, deltaT --> Timestemp
	xLen = XDIM / vproc[0]; // for extra rows from the other processors
	yLen = YDIM / vproc[1];
	zLen = ZDIM / vproc[2];
	
	setNeighboursAndCommunicationConsts();

	time(&startTime);
	memoryAllocation(xLen + 2, yLen + 2, zLen + 2);
	time(&endTime);
	runTime = difftime(endTime,startTime);	
	//printf("\n processor # %d memoryAllocation %e \n",sid, runTime);
	// set all the materials and sources
	time(&startTime);
	defMaterial();// define all the material types
	time(&endTime);	
	runTime = difftime(endTime,startTime);
	//printf("\n processor # %d defMaterial%e \n",sid, runTime);
			
	time(&startTime);
	defSources();// define all the sources 
	time(&endTime);
	runTime = difftime(endTime,startTime);
	//printf("\n processor # %d defSources %e \n",sid, runTime);
	// No need to initialize sources as update sources are called before each step (including the first step)
	
	// Ca and Da are dependant on the material properties and hence geometry

	time(&startTime);
	defGeometry();
	time(&endTime);
	runTime = difftime(endTime,startTime);
	printf("\n processor # %d defGeometry %e \n",sid, runTime);
	
	// for PML
	if (PML_LAYERS!= 0)
		allocatePMLFields();
	allocateDispFields();				
} // end init
