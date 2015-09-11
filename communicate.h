void createSendBuffer(int xStart, int xEnd, int yStart, int yEnd, int zStart, int zEnd, double *buffer, int fieldIndex1, int fieldIndex2);
void readRecvBuffer(int xStart, int xEnd, int yStart, int yEnd, int zStart, int zEnd, double *buffer, int fieldIndex1, int fieldIndex2);

void communicateHfield(int n);
void communicateHfieldNS(int n);
void communicateHfieldEW(int n);
void communicateHfieldTB(int n);


void communicateEfield(int n);
void communicateEfieldEW(int n);
void communicateEfieldNS(int n);
void communicateEfieldTB(int n);

