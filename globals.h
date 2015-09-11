///////////////////
//    Globals     /
// Hence init to 0/  
///////////////////
// For the Ex, Ey and Ez etc and material index

fieldMaterialStruct ***fieldMat;

sourceStruct source[MAXMATERIALS];
int listOfSources[MAXMATERIALS]; // Stores the material number of all the sources


materialStruct material[MAXMATERIALS];

double deltaT;

double sigmaNot;

int coordinates[3];
int xLen, yLen, zLen;// number of points of the present computer only, add 2 to get the points from the other rows
int px, py, pz; // processors coordinates
int sid, vproc[3];

// For communication
double *sendBuffer1;
double *recvBuffer1;

double *sendBuffer2;
double *recvBuffer2;
MPI_Comm grid_comm;

// neighbours
int nebW,nebE,nebS,nebN,nebB,nebT;	
int copyCountEW, copyCountNS, copyCountTB;

pmlConstStruct cn[PML_LAYERS];
pmlConstStruct cOnes[PML_LAYERS]; 

pmlFieldStruct *pmlFields;


lorentzParam lorentz[10];

int dispFieldIndex;
dispFieldStruct *dispField;
int totalNumberOfMaterials;
