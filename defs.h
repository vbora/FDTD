///////////////////
//    Globals     /
// Hence init to 0/  
///////////////////
// For the Ex, Ey and Ez etc and material index

typedef struct{
	double field[6];	
	int materialIndex;
	int pmlIndex;
	int dispFieldIndex; //Index in dispField (Stored only E fields)
} fieldMaterialStruct;


typedef struct{
	double dx;
	double dy;
	double dz;
	double bx;
	double by;
	double bz;
} pmlFieldStruct;


typedef struct{
	double muR;
	double epsR;
	double sigma;
	double Ca;
	double Cb;
	double Da;
	double Db; 
	double Jsx;
	double Jsy;
	double Jsz;
	double Msx;
	double Msy;
	double Msz;
        int    numberOfPoles; // Total number of poles
        int    lorentzIndex;  // The index of the first pole in lorentz matrix
	double dispC1;
	double dispC2;
	double dispC3;
} materialStruct;

typedef struct{
	int sourceType; // Like continuos, pulse, gaussian
	double frequency;
	double deltaFrequency;
	double phase;
	int delay; // in time steps
	int polarization;
	double peakValue;
	int location; 	// This is the material index
		 	// So this will traverse through all the grid points!	
} sourceStruct;

typedef struct{
	double c1;
	double c2;
	double c3;
	double c4;
	double c5;
	double c6;
} pmlConstStruct;

typedef struct{
	double deltaEpsilon;
	double omegaPlasmon;
	double deltaPlasmon;
	
	double alpha;
	double zeta;
	double gamma;
	
} lorentzParam;
/*
typedef struct{
	double exMinus1;
	double eyMinus1;
	double ezMinus1;
	int dispCurrentIndex; // index of the first pole in dispCurrent
} dispFieldStruct;

typedef struct{
	double jdx;
	double jdxMinus1;
	double jdy;
	double jdyMinus1;
	double jdz;
	double jdzMinus1;
} dispCurrentStruct;
*/

typedef struct{
	double exMinus1; // e's used only for the first pole, so wasting space because of allocation deadlock
	double eyMinus1;
	double ezMinus1;

	double jdx;
	double jdxMinus1;
	double jdy;
	double jdyMinus1;
	double jdz;
	double jdzMinus1;
} dispFieldStruct;
