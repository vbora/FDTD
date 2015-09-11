///////////////////
//    Globals     /
// Hence init to 0/  
///////////////////
// For the Ex, Ey and Ez

extern fieldMaterialStruct ***fieldMat;

extern sourceStruct source[MAXMATERIALS];
extern int listOfSources[MAXMATERIALS]; // Stores the material number of all the sources


extern materialStruct material[MAXMATERIALS];

extern double deltaT;

extern double sigmaNot;

extern int coordinates[3];
extern int xLen, yLen, zLen;// number of points of the present computer only, add 2 to get the points from the other rows
extern int px, py, pz; // processors coordinates
extern int sid, vproc[3];

// For communication
extern double *sendBuffer1;
extern double *recvBuffer1;

extern double *sendBuffer2;
extern double *recvBuffer2;
extern MPI_Comm grid_comm;	
// neighbours
extern int nebW,nebE,nebS,nebN,nebB,nebT;
extern int copyCountEW, copyCountNS, copyCountTB;

extern pmlConstStruct cn[PML_LAYERS];
extern pmlConstStruct cOnes[PML_LAYERS];
extern pmlFieldStruct *pmlFields;

extern lorentzParam lorentz[10];

extern int dispFieldIndex;
dispFieldStruct *dispField;
extern int totalNumberOfMaterials;

