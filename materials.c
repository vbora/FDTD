#include"common.h"
#include"defs.h"
#include"externGlobal.h"
#include"materials.h"

void setOrdinaryMaterial(int m){
		CA(m) = (1 - SIGMA(m)*deltaT/(2*EPSR(m)*EPSNOT )) / ( 1 + SIGMA(m)*deltaT /(2*EPSR(m)*EPSNOT));   
		CB(m) = (deltaT/(EPSR(m)*EPSNOT*DELTA)) / ( 1 + SIGMA(m) * deltaT /(2*EPSR(m)*EPSNOT));   
		DA(m) = (1 - SIGMA(m)*deltaT/(2*MUR(m)*MUNOT)) / ( 1 + SIGMA(m) * deltaT /(2*MUR(m)*MUNOT));   
		DB(m) = (deltaT/(MUR(m)*MUNOT*DELTA)) / ( 1 + SIGMA(m)*deltaT /(2*MUR(m)*MUNOT));	
} // end setOrdinaryMaterial
