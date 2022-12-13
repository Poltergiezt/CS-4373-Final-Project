#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>


double determinantOfMatrix(double* matrix, int matrixSize);
double logDeterminantOfMatrix(double* matrix, int matrixSize);
double timeDiff(struct timeval* start, struct timeval* end);


int main(int argc, char* argv[])
{
	struct timeval programStart;
	gettimeofday(&programStart, NULL);
	if(argc < 2)
	{
		fprintf(stderr, "Usage: %s matrix_size\n", argv[0]);
		exit(1);
	}
	int matrixSize = strtol(argv[1], NULL, 10);
	char filePath[50];
	sprintf(filePath, "../DET_DATA/");
	char matrixSizeString[10];
	sprintf(matrixSizeString, "%d", matrixSize);
	char paddedMatrixSizeString[10];
	char fileName[50]; // Construct file name separately because problem requirements demand it.
	sprintf(fileName, "m");
	if(matrixSize > 999)
	{
		sprintf(paddedMatrixSizeString, "%s", matrixSizeString);
	}
	else if(matrixSize > 99)
	{
		sprintf(paddedMatrixSizeString, "0");
		strcat(paddedMatrixSizeString, matrixSizeString);
	}
	else
	{
		sprintf(paddedMatrixSizeString, "00");
		strcat(paddedMatrixSizeString, matrixSizeString);
	}
	strcat(fileName, paddedMatrixSizeString);
	strcat(fileName, "x");
	strcat(fileName, paddedMatrixSizeString);
	strcat(fileName, ".bin");
	strcat(filePath, fileName);
	FILE* dataFile = fopen(filePath, "rb");
	if(dataFile == NULL)
	{
		fprintf(stderr, "File %s does not exist.\n", filePath);
		exit(1);
	}
	printf("(1) %s\n", fileName);
	printf("(2) %d\n", matrixSize);
	// Read elements
	double* matrix = calloc(matrixSize * matrixSize, sizeof(double));
	for(int i = 0; i < matrixSize; i++)
	{
		for(int j = 0; j < matrixSize; j++)
		{
			fread(&matrix[i * matrixSize + j], sizeof(double), 1, dataFile);
		}
	}
	fclose(dataFile);

	struct timeval workStart;
	gettimeofday(&workStart, NULL);
	// Calculate determinant
	double det = determinantOfMatrix(matrix, matrixSize);
	double logDet = logDeterminantOfMatrix(matrix, matrixSize);
	free(matrix);
	struct timeval workEnd;
	gettimeofday(&workEnd, NULL);
	struct timeval programEnd;
	gettimeofday(&programEnd, NULL);

	printf("(3) %f\n", det);
	printf("(4) %f\n", logDet);
	printf("%f\n", timeDiff(&workStart, &workEnd));
	printf("%f\n", timeDiff(&programStart, &programEnd));

}


double determinantOfMatrix(double* matrix, int matrixSize)
{
	int i, j, k;
	double det = 1;
	double ratio;

	/* Applying Gauss Elimination */
	for(i = 0; i < matrixSize; i++)
	{
		if(matrix[i * matrixSize + i] == 0.0)
		{
			printf("Mathematical Error!");
			exit(0);
		}
		for(j = i + 1; j < matrixSize; j++)
		{
			ratio = matrix[j * matrixSize + i] / matrix[i * matrixSize + i];

			for(k = 0; k < matrixSize; k++)
			{
				matrix[j * matrixSize + k] = matrix[j * matrixSize + k] - ratio * matrix[i * matrixSize + k];
			}
		}
	}

	/* Finding determinant by multiplying
	 elements in principal diagonal elements */
	for(i = 0; i < matrixSize; i++)
	{
		det = det * matrix[i * matrixSize + i];
	}

	return det;
}


double logDeterminantOfMatrix(double* matrix, int matrixSize)
{
	int i, j, k;
	double logDet = 0;
	double ratio;

	/* Applying Gauss Elimination */
	for(i = 0; i < matrixSize; i++)
	{
		if(matrix[i * matrixSize + i] == 0.0)
		{
			printf("Mathematical Error!");
			exit(0);
		}
		for(j = i + 1; j < matrixSize; j++)
		{
			ratio = matrix[j * matrixSize + i] / matrix[i * matrixSize + i];

			for(k = 0; k < matrixSize; k++)
			{
				matrix[j * matrixSize + k] = matrix[j * matrixSize + k] - ratio * matrix[i * matrixSize + k];
			}
		}
	}

	/* Finding determinant by multiplying
	 elements in principal diagonal elements */
	for(i = 0; i < matrixSize; i++)
	{
		logDet = logDet + log10(fabs(matrix[i * matrixSize + i]));
	}

	return logDet;
}


double timeDiff(struct timeval* start, struct timeval* end)
{
	return ((double) end->tv_sec - (double) start->tv_sec) + 1e-6 * (end->tv_usec - start->tv_usec);
}