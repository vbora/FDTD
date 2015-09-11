#include"common.h"
#include"defs.h"
#include"externGlobal.h"
#include"fileUtil.h"

void writeMatFile(char baseName[20]){
        int ix,iy, iz, m;
        char delim, fileName[50];
        FILE *outFile;
                                                                                
        delim = ',';
                                                                                
        sprintf(fileName,"%s%d",baseName,sid);
        //strcpy(fileName, baseName);
        //strcat(fileName,sNum);
        outFile = fopen(fileName,"w");
	m = MATERIALINDEX(1,1,1);
        //fprintf(outFile, "%e", EPSR(m));
       	for( iz=1; iz < zLen + 1; iz++)
       	for( iy=1; iy < yLen + 1; iy++)
	for( ix=1; ix < xLen + 1; ix++){
		m = MATERIALINDEX(ix,iy,iz);
       		 fprintf(outFile, "%e%c",EPSR(m), delim);
        }
        fclose(outFile);
} // end writeFile

void writeFile(int fieldComponent, char baseName[20], int stepNumber){
        int ix,iy, iz;
        char delim, fileName[50];
        FILE *outFile;
                                                                                
        delim = ',';
                                                                                
        sprintf(fileName,"%d%s%d", stepNumber, baseName,sid);
        //strcpy(fileName, baseName);
        //strcat(fileName,sNum);
        outFile = fopen(fileName,"w");
       	for( iz=1; iz < zLen + 1; iz++)
        for( iy=1; iy < yLen + 1; iy++)
	for( ix=1; ix < xLen + 1; ix++){
        	                fprintf(outFile, "%e%c", FIELD(ix,iy,iz,fieldComponent), delim);
        }
	
        fclose(outFile);
} // end writeFile

void printStep(int n){
        //if ((n > 400) && (n % 20 == 0))
        char filePrefix[100];
        if (n % 100 == 0 && n != 0){
	writeFile(HZOFT, "HZ", n);
	writeFile(EYOFT, "EY", n);
/*                sprintf(filePrefix,"%d",sid);
                strcat(filePrefix,"Hz");
                writeFile(Hz,filePrefix,n);
                                                                                
                sprintf(filePrefix,"%d",sid);
                strcat(filePrefix,"Ex");
                writeFile(Ex,filePrefix,n);
                                                                                
                sprintf(filePrefix,"%d",sid);
                strcat(filePrefix,"Ey");
                writeFile(Ey,filePrefix,n);
*/
        }
                                                                                
}// end printstep


