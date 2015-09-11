#include<stdio.h>
#include<time.h>
#include"mpi.h"
#include"input.h"

//####
//Constants
//Moved to input.h
/*
#define XDIM 1000
#define YDIM 500
#define ZDIM 200
#define TOTALTIMESTEPS 5000
#define PML_LAYERS 10
*/ 
//#define MUNOT 1.25663706e-6
//#define EPSNOT  8.85418782e-12
//#define c 299792458
// The constants have to be normalized since c is reduced to 1
// Taking Dr. Prata's c = 299792458
#define c 1 
#define MUNOT 376.73031303129346
#define EPSNOT  0.002654418730151
#define MAXMATERIALS 10
#define PI 3.14159265358979323846


// Courant factor
#define S 0.5
// Assuming deltaX = deltaY = deltaZ = DELTA, more about this soon!TODO 
#define DELTA 1

#define MUR(m)  material[m].muR
#define EPSR(m) material[m].epsR
#define SIGMA(m) material[m].sigma
#define CA(m) 	material[m].Ca
#define CB(m) 	material[m].Cb
#define DA(m) 	material[m].Da
#define DB(m)	material[m].Db
#define JSX(m)	material[m].Jsx 
#define JSY(m)	material[m].Jsy 
#define JSZ(m)	material[m].Jsz 
#define MSX(m)	material[m].Msx 
#define MSY(m)	material[m].Msy 
#define MSZ(m)	material[m].Msz 
#define NPOLES(m) material[m].numberOfPoles
#define LORN_INDEX(m) material[m].lorentzIndex
#define DISPC1(m) material[m].dispC1
#define DISPC2(m) material[m].dispC2
#define DISPC3(m) material[m].dispC3

// PML material data structure

#define EXOFT 0
#define EYOFT 1
#define EZOFT 2
#define HXOFT 3
#define HYOFT 4
#define HZOFT 5

#define MATERIALINDEX(i,j,k) fieldMat[i][j][k].materialIndex 
#define EX(i,j,k) fieldMat[i][j][k].field[EXOFT] 
#define EY(i,j,k) fieldMat[i][j][k].field[EYOFT] 
#define EZ(i,j,k) fieldMat[i][j][k].field[EZOFT] 
#define HX(i,j,k) fieldMat[i][j][k].field[HXOFT] 
#define HY(i,j,k) fieldMat[i][j][k].field[HYOFT] 
#define HZ(i,j,k) fieldMat[i][j][k].field[HZOFT] 
#define PMLINDEX(i,j,k) fieldMat[i][j][k].pmlIndex 
#define DISPFIELDINDEX(i,j,k) fieldMat[i][j][k].dispFieldIndex 

// Alternate way of accessing the fields, this is useful for printing and stuff like that
#define FIELD(i,j,k,type) fieldMat[i][j][k].field[type]
 
// sourceStruct elements
#define SOURCETYPE(m) source[m].sourceType
#define FREQUENCY(m) source[m].frequency
#define DELTA_FREQUENCY(m) source[m].frequency
#define PHASE(m) source[m].phase
#define DELAY(m) source[m].delay
#define POLARIZATION(m) source[m].polarization
#define PEAKVALUE(m) source[m].peakValue
#define LOCATION(m) source[m].location


// Source constants also used for writing files

// polarization
#define EXFIELD 1 // No default, will give error if not stated
#define EYFIELD 2
#define EZFIELD 3

#define HXFIELD 4
#define HYFIELD 5
#define HZFIELD 6

// Source type
#define PULSE 0 // Default
#define SINE 1 
#define GAUSSIAN 2 

// Geometry features
#define X 0 
#define Y 1
#define Z 2

// Edges names
/*
        North
         __
West    |  | East
        |__|
  (0,0) South

*/
#define NORTH  0
#define WEST   1
#define SOUTH  2
#define EAST   3
#define TOP    4
#define BOTTOM 5

// for split problem
#define MAX_FACT 10000

////////////////////////
// PML layers

// Grading of the PML
// Choosing polynomial grading
#define REFL0 1e-6

// value of  3<=EXP_M<=4
#define EXP_M 3
#define ETTA  376.7303131959
#define KAPPA_MAX 100

#define DX(m) pmlFields[m].dx
#define DY(m) pmlFields[m].dy
#define DZ(m) pmlFields[m].dz
#define BX(m) pmlFields[m].bx
#define BY(m) pmlFields[m].by
#define BZ(m) pmlFields[m].bz


// Dispersive Material
#define DELTAEPS(p) lorentz[p].deltaEpsilon
#define OMEGAP(p)   lorentz[p].omegaPlasmon
#define DELTAP(p)   lorentz[p].deltaPlasmon

#define ALPHA(p)    lorentz[p].alpha
#define ZETA(p)     lorentz[p].zeta
#define GAMMA(p)    lorentz[p].gamma


#define JDX(d)       dispField[d].jdx
#define JDXMINUS1(d) dispField[d].jdxMinus1
#define JDY(d)       dispField[d].jdy
#define JDYMINUS1(d) dispField[d].jdyMinus1
#define JDZ(d)       dispField[d].jdz
#define JDZMINUS1(d) dispField[d].jdzMinus1

#define EXMINUS1(d)  dispField[d].exMinus1
#define EYMINUS1(d)  dispField[d].eyMinus1
#define EZMINUS1(d)  dispField[d].ezMinus1
