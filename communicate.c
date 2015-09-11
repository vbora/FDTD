#include"common.h"
#include"defs.h"
#include"externGlobal.h"
#include"util.h"
#include"communicate.h"

void communicateEfield(int n){
	//printf("sid %d Neighbours nebW %d, nebE %d, nebS %d, nebN %d, nebT%d, nebB %d\n",sid, nebW, nebE, nebS, nebN, nebT, nebB);
	communicateEfieldEW(n);	
	//printf("communicateEfieldEW processor # %d\n ",sid );	
	communicateEfieldNS(n);	
	//printf("communicateEfieldNS processor # %d\n ",sid );	
	communicateEfieldTB(n);	
	//printf("communicateEfieldTB processor # %d\n ",sid );	
	//printf("sid %d  after E communicate Neighbours nebW %d, nebE %d, nebS %d, nebN %d, nebT%d, nebB %d\n",sid, nebW, nebE, nebS, nebN, nebT, nebB);
	MPI_Barrier(MPI_COMM_WORLD);
	//printf("Finished E communication processor %d\n", sid);

} // end communicateEfield 

void communicateHfield(int n){
	communicateHfieldEW(n);
	//printf("communicateHfieldEW processor # %d\n ",sid );	
	communicateHfieldNS(n);	
	//printf("communicateHfieldNS processor # %d\n ",sid );	
	communicateHfieldTB(n);	
	//printf("communicateHfieldTB processor # %d\n ",sid );	
	//printf("sid %d  after H communicate Neighbours nebW %d, nebE %d, nebS %d, nebN %d, nebT%d, nebB %d\n",sid, nebW, nebE, nebS, nebN, nebT, nebB);
	MPI_Barrier(grid_comm);
	//printf("Finished H communication processor %d\n", sid);

} // end communicateHfield 

void createSendBuffer(int xStart, int xEnd, int yStart, int yEnd, int zStart, int zEnd, double *buffer, int fieldIndex1, int fieldIndex2){
	int index = 0;
	int i,j,k;	
	for(i=xStart; i<xEnd+1; i++)
	for(j=yStart; j<yEnd+1; j++)
	for(k=zStart; k<zEnd+1; k++){
		buffer[index]   = FIELD(i,j,k,fieldIndex1); 	
		buffer[index+1] = FIELD(i,j,k,fieldIndex2); 	
		index += 2; 
	}
} // end function createSendBuffer

void readRecvBuffer(int xStart, int xEnd, int yStart, int yEnd, int zStart, int zEnd, double *buffer, int fieldIndex1, int fieldIndex2){
	int index = 0;
	int i,j,k;	
	for(i=xStart; i<xEnd+1; i++)
	for(j=yStart; j<yEnd+1; j++)
	for(k=zStart; k<zEnd+1; k++){
		FIELD(i,j,k,fieldIndex1) = buffer[index]; 	
		FIELD(i,j,k,fieldIndex2) = buffer[index+1]; 	
		index += 2; 
	}

} // end function readRecvBuffer 




void communicateHfieldEW(int n){

	// since EW, no normal component ie Hx
	int i,j,k,index;
	MPI_Request recvRequestE, recvRequestW;	// The E/W alwayes points to the neighbour, be careful with the tags	

	MPI_Status recvStatusE, recvStatusW;

	double *sendBufferW, *sendBufferE, *recvBufferW, *recvBufferE;

	sendBufferE = sendBuffer1;
	sendBufferW = sendBuffer2;

	recvBufferE = recvBuffer1;
	recvBufferW = recvBuffer2;
	//printf("Processor %d nebW %d, nebE %d, nebN %d, nebS %d, nebT %d, nebB %d\n",sid, nebW, nebE, nebN, nebS, nebT, nebB);
	// async recv
	if (nebW != -1)		
		MPI_Irecv( recvBufferW, copyCountEW , MPI_DOUBLE, nebW,10, grid_comm, &recvRequestW );
	if (nebE != -1)		
		MPI_Irecv( recvBufferE, copyCountEW , MPI_DOUBLE, nebE,20, grid_comm, &recvRequestE );

		//send
	
	if (nebW != -1){		
		createSendBuffer(1, 1, 1, yLen, 1, zLen, sendBufferW, HYOFT,HZOFT);
		MPI_Send( sendBufferW, copyCountEW, MPI_DOUBLE, nebW, 20, grid_comm);
	}

	if (nebE != -1){		
		createSendBuffer(xLen, xLen, 1, yLen, 1, zLen, sendBufferE, HYOFT,HZOFT);
		MPI_Send( sendBufferE, copyCountEW, MPI_DOUBLE, nebE, 10, grid_comm);
	}
	
	if (nebW != -1){
		MPI_Wait(&recvRequestW,&recvStatusW);
		readRecvBuffer(0, 0, 1, yLen, 1, zLen, recvBufferW, HYOFT, HZOFT);
	}
		
	if (nebE != -1){
		MPI_Wait(&recvRequestE,&recvStatusE);
		readRecvBuffer(xLen+1, xLen+1, 1, yLen, 1, zLen, recvBufferE, HYOFT, HZOFT);
	}		
}


void communicateHfieldNS(int n){

	// since NS, no normal component ie Hy
	int i,j,k,index;
	MPI_Request recvRequestS, recvRequestN;	// The E/W alwayes points to the neighbour, be careful with the tags	

	MPI_Status recvStatusS, recvStatusN;

	double *sendBufferS, *sendBufferN, *recvBufferS, *recvBufferN;

	sendBufferS = sendBuffer1;
	sendBufferN = sendBuffer2;

	recvBufferS = recvBuffer1;
	recvBufferN = recvBuffer2;

	// async recv
	if (nebS != -1)		
		MPI_Irecv( recvBufferS, copyCountNS , MPI_DOUBLE, nebS,10, grid_comm, &recvRequestS );
	if (nebN != -1)		
		MPI_Irecv( recvBufferN, copyCountNS , MPI_DOUBLE, nebN,20, grid_comm, &recvRequestN );

		//send
	
	if (nebS != -1){		
		createSendBuffer(1, xLen, 1, 1, 1, zLen, sendBufferS, HXOFT,HZOFT);
		MPI_Send( sendBufferS, copyCountNS, MPI_DOUBLE, nebS, 20, grid_comm);
	}

	if (nebN != -1){		
		createSendBuffer(1, xLen, yLen, yLen, 1, zLen, sendBufferN, HXOFT,HZOFT);
		MPI_Send( sendBufferN, copyCountNS, MPI_DOUBLE, nebN, 10, grid_comm);
	}
	
	if (nebS != -1){
		MPI_Wait(&recvRequestS,&recvStatusS);
		readRecvBuffer(1, xLen, 0, 0, 1, zLen, recvBufferS, HXOFT, HZOFT);
	}
		
	if (nebN != -1){
		MPI_Wait(&recvRequestN,&recvStatusN);
		readRecvBuffer(1, xLen, yLen+1, yLen+1, 1, zLen, recvBufferN, HXOFT, HZOFT);
	}		
}// end communicateHfieldNS



void communicateHfieldTB(int n){

	// since TB, no normal component ie Hz
	int i,j,k,index;
	MPI_Request recvRequestB, recvRequestT;	// The E/W alwayes points to the neighbour, be careful with the tags	

	MPI_Status recvStatusB, recvStatusT;

	double *sendBufferB, *sendBufferT, *recvBufferB, *recvBufferT;

	sendBufferB = sendBuffer1;
	sendBufferT = sendBuffer2;

	recvBufferB = recvBuffer1;
	recvBufferT = recvBuffer2;

	// async recv
	if (nebB != -1)		
		MPI_Irecv( recvBufferB, copyCountTB , MPI_DOUBLE, nebB,10, grid_comm, &recvRequestB );
	if (nebT != -1)		
		MPI_Irecv( recvBufferT, copyCountTB , MPI_DOUBLE, nebT,20, grid_comm, &recvRequestT );

	//send
	
	if (nebB != -1){		
		createSendBuffer(1, xLen, 1, yLen, 1, 1, sendBufferB, HXOFT,HYOFT);
		MPI_Send( sendBufferB, copyCountTB, MPI_DOUBLE, nebB, 20, grid_comm);
	}

	if (nebT != -1){		
		createSendBuffer(1, xLen, 1, yLen, zLen, zLen, sendBufferT, HXOFT,HYOFT);
		MPI_Send( sendBufferT, copyCountTB, MPI_DOUBLE, nebT, 10, grid_comm);
	}
	// fill buffer	
	if (nebB != -1){
		MPI_Wait(&recvRequestB,&recvStatusB);
		readRecvBuffer(1, xLen, 1, yLen, 0, 0, recvBufferB, HXOFT, HYOFT);
	}
		
	if (nebT != -1){
		MPI_Wait(&recvRequestT,&recvStatusT);
		readRecvBuffer(1, xLen, 1, yLen, zLen+1, zLen+1, recvBufferT, HXOFT, HYOFT);
	}		
} // end communicateHfieldTB


//////////////////////////////////E fields

void communicateEfieldEW(int n){

	// since EW, no normal component ie Ex
	int i,j,k,index;
	MPI_Request recvRequestE, recvRequestW;	// The E/W alwayes points to the neighbour, be careful with the tags	

	MPI_Status recvStatusE, recvStatusW;

	double *sendBufferW, *sendBufferE, *recvBufferW, *recvBufferE;

	sendBufferE = sendBuffer1;
	sendBufferW = sendBuffer2;

	recvBufferE = recvBuffer1;
	recvBufferW = recvBuffer2;

	// async recv
	if (nebW != -1)		
		MPI_Irecv( recvBufferW, copyCountEW , MPI_DOUBLE, nebW,10, grid_comm, &recvRequestW );
	if (nebE != -1)		
		MPI_Irecv( recvBufferE, copyCountEW , MPI_DOUBLE, nebE,20, grid_comm, &recvRequestE );

		//send
	
	if (nebW != -1){		
		createSendBuffer(1, 1, 1, yLen, 1, zLen, sendBufferW, EYOFT,EZOFT);
		MPI_Send( sendBufferW, copyCountEW, MPI_DOUBLE, nebW, 20, grid_comm);
	}

	if (nebE != -1){		
		createSendBuffer(xLen, xLen, 1, yLen, 1, zLen, sendBufferE, EYOFT,EZOFT);
		MPI_Send( sendBufferE, copyCountEW, MPI_DOUBLE, nebE, 10, grid_comm);
	}
	
	if (nebW != -1){
		MPI_Wait(&recvRequestW,&recvStatusW);
		readRecvBuffer(0, 0, 1, yLen, 1, zLen, recvBufferW, EYOFT, EZOFT);
	}
		
	if (nebE != -1){
		MPI_Wait(&recvRequestE,&recvStatusE);
		readRecvBuffer(xLen+1, xLen+1, 1, yLen, 1, zLen, recvBufferE, EYOFT, EZOFT);
	}		
}


void communicateEfieldNS(int n){

	// since NS, no normal component ie Ey
	int i,j,k,index;
	MPI_Request recvRequestS, recvRequestN;	// The E/W alwayes points to the neighbour, be careful with the tags	

	MPI_Status recvStatusS, recvStatusN;

	double *sendBufferS, *sendBufferN, *recvBufferS, *recvBufferN;

	sendBufferS = sendBuffer1;
	sendBufferN = sendBuffer2;

	recvBufferS = recvBuffer1;
	recvBufferN = recvBuffer2;

	// async recv
	if (nebS != -1)		
		MPI_Irecv( recvBufferS, copyCountNS , MPI_DOUBLE, nebS,10, grid_comm, &recvRequestS );
	if (nebN != -1)		
		MPI_Irecv( recvBufferN, copyCountNS , MPI_DOUBLE, nebN,20, grid_comm, &recvRequestN );

		//send
	
	if (nebS != -1){		
		createSendBuffer(1, xLen, 1, 1, 1, zLen, sendBufferS, EXOFT,EZOFT);
		MPI_Send( sendBufferS, copyCountNS, MPI_DOUBLE, nebS, 20, grid_comm);
	}

	if (nebN != -1){		
		createSendBuffer(1, xLen, yLen, yLen, 1, zLen, sendBufferN, EXOFT,EZOFT);
		MPI_Send( sendBufferN, copyCountNS, MPI_DOUBLE, nebN, 10, grid_comm);
	}
	
	if (nebS != -1){
		MPI_Wait(&recvRequestS,&recvStatusS);
		readRecvBuffer(1, xLen, 0, 0, 1, zLen, recvBufferS, EXOFT, EZOFT);
	}
		
	if (nebN != -1){
		MPI_Wait(&recvRequestN,&recvStatusN);
		readRecvBuffer(1, xLen, yLen+1, yLen+1, 1, zLen, recvBufferN, EXOFT, EZOFT);
	}		
}// end communicateEfieldNS



void communicateEfieldTB(int n){

	// since TB, no normal component ie Ez
	int i,j,k,index;
	MPI_Request recvRequestB, recvRequestT;	// The E/W alwayes points to the neighbour, be careful with the tags	

	MPI_Status recvStatusB, recvStatusT;

	double *sendBufferB, *sendBufferT, *recvBufferB, *recvBufferT;

	sendBufferB = sendBuffer1;
	sendBufferT = sendBuffer2;

	recvBufferB = recvBuffer1;
	recvBufferT = recvBuffer2;

	// async recv
	if (nebB != -1)		
		MPI_Irecv( recvBufferB, copyCountTB , MPI_DOUBLE, nebB,10, grid_comm, &recvRequestB );
	if (nebT != -1)		
		MPI_Irecv( recvBufferT, copyCountTB , MPI_DOUBLE, nebT,20, grid_comm, &recvRequestT );

	//send
	
	if (nebB != -1){		
		createSendBuffer(1, xLen, 1, yLen, 1, 1, sendBufferB, EXOFT,EYOFT);
		MPI_Send( sendBufferB, copyCountTB, MPI_DOUBLE, nebB, 20, grid_comm);
	}

	if (nebT != -1){		
		createSendBuffer(1, xLen, 1, yLen, zLen, zLen, sendBufferT, EXOFT,EYOFT);
		MPI_Send( sendBufferT, copyCountTB, MPI_DOUBLE, nebT, 10, grid_comm);
	}
	// fill buffer	
	if (nebB != -1){
		MPI_Wait(&recvRequestB,&recvStatusB);
		readRecvBuffer(1, xLen, 1, yLen, 0, 0, recvBufferB, EXOFT, EYOFT);
	}
		
	if (nebT != -1){
		MPI_Wait(&recvRequestT,&recvStatusT);
		readRecvBuffer(1, xLen, 1, yLen, zLen+1, zLen+1, recvBufferT, EXOFT, EYOFT);
	}		
} // end communicateEfieldTB

