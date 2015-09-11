#include"common.h"
#include"defs.h"
#include"externGlobal.h"

void updateSource(int n){
	int numberOfSources;
	int i,m;
	double scalarValue;
	i = 0;
	while(listOfSources[i] != -1){
		m = listOfSources[i]; // also equal to LOCATION_M
		//printf(" sid %d m %d\n",sid,m);
	 	switch(SOURCETYPE(m)){
			case PULSE: {
				if (n - DELAY(m) ==0){
					if (POLARIZATION(m) == EXFIELD)
						JSX(m) = PEAKVALUE(m);
					else if(POLARIZATION(m) == EYFIELD)
						JSY(m) = PEAKVALUE(m);
					else if(POLARIZATION(m) == EZFIELD)
						JSZ(m) = PEAKVALUE(m);
					else if(POLARIZATION(m) == HXFIELD)
						MSX(m) = PEAKVALUE(m);
					else if(POLARIZATION(m) == HYFIELD)
						MSY(m) = PEAKVALUE(m);
					else if(POLARIZATION(m) == HZFIELD)
						MSZ(m) = PEAKVALUE(m);
				}	
				else{
					JSX(m) = 0;		
					JSY(m) = 0;		
					JSZ(m) = 0;		
					MSX(m) = 0;		
					MSY(m) = 0;		
					MSZ(m) = 0;		
				}
				break;
			} // end PULSE	
			case SINE:{
							//calculate scalar value and then add polarization
					scalarValue = PEAKVALUE(m) * sin(2*PI* FREQUENCY(m) * n + PHASE(m)/(FREQUENCY(m)*2*PI));		
					//printf("\n processors # %d Source value on time step %d -->  %e\n",sid,n, scalarValue); 	
					if (POLARIZATION(m) == EXFIELD)
						JSX(m) = scalarValue;
					else if(POLARIZATION(m) == EYFIELD)
						JSY(m) = scalarValue;
					else if(POLARIZATION(m) == EZFIELD)
						JSZ(m) = scalarValue;
					else if(POLARIZATION(m) == HXFIELD)
						MSX(m) = scalarValue;
					else if(POLARIZATION(m) == HYFIELD)
						MSY(m) = scalarValue;
					else if(POLARIZATION(m) == HZFIELD)
						MSZ(m) = scalarValue;
					break;
			}// end SINE
	}// end case
	i++;	
	}// end while	
}// end update source

