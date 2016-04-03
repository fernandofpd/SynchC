/* Fernando Perez-Diaz, February 2016           */

int readInputFile(char *inputFilename, double **matrix, int numRows, int numCols, double minVal, double maxVal);
void inputOptions(int argc, char *argv[], int *outConn, int *outOrdPar, int *outInterspike, char *ftsync, char  *fphase, char *ftau, char *fpos, char *fvel);
