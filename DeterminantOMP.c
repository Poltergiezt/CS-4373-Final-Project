#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>


double determinantOfMatrix(double* matrix, int matrixSize);
double logDeterminantOfMatrix(double* matrix, int matrixSize);


int g_threadCount;


int main(int argc, char* argv[])
{
	double programStart = omp_get_wtime();
	if(argc < 3)
	{
		fprintf(stderr, "Usage: %s matrix_size thread_count\n", argv[0]);
		exit(1);
	}
	int matrixSize = strtol(argv[1], NULL, 10);
	g_threadCount = strtol(argv[2], NULL, 10);
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

	double workStart = omp_get_wtime();
	// Calculate determinant
	double det = determinantOfMatrix(matrix, matrixSize);
	double logDet = logDeterminantOfMatrix(matrix, matrixSize);
	free(matrix);
	double workEnd = omp_get_wtime();
	double programEnd = omp_get_wtime();

	printf("(3) %.6e\n", det);
	printf("(4) %.6e\n", logDet);
	printf("%f\n", workEnd - workStart);
	printf("%f\n", programEnd - programStart);
}


double determinantOfMatrix(double* matrix, int matrixSize)
{
	int i, j;
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
		#pragma omp parallel for num_threads(g_threadCount) default(none) private(ratio) shared(matrix, matrixSize, i)
		for(j = i + 1; j < matrixSize; j++)
		{
			ratio = matrix[j * matrixSize + i] / matrix[i * matrixSize + i];

			for(int k = 0; k < matrixSize; k++)
			{
				matrix[j * matrixSize + k] = matrix[j * matrixSize + k] - ratio * matrix[i * matrixSize + k];
			}
		}
	}

	/* Finding determinant by multiplying
	 elements in principal diagonal elements */
	#pragma omp parallel for num_threads(g_threadCount) default(none) reduction(*: det) shared(matrix, matrixSize)
	for(i = 0; i < matrixSize; i++)
	{
		det *= matrix[i * matrixSize + i];
	}

	return det;
}


double logDeterminantOfMatrix(double* matrix, int matrixSize)
{
	int i, j;
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
		#pragma omp parallel for num_threads(g_threadCount) private(ratio) shared(matrix, matrixSize, i) default(none)
		for(j = i + 1; j < matrixSize; j++)
		{
			ratio = matrix[j * matrixSize + i] / matrix[i * matrixSize + i];

			for(int k = 0; k < matrixSize; k++)
			{
				matrix[j * matrixSize + k] = matrix[j * matrixSize + k] - ratio * matrix[i * matrixSize + k];
			}
		}
	}

	/* Finding determinant by multiplying
	 elements in principal diagonal elements */
	#pragma omp parallel for num_threads(g_threadCount) reduction(+: logDet) shared(matrix, matrixSize) default(none)
	for(i = 0; i < matrixSize; i++)
	{
		logDet += log10(fabs(matrix[i * matrixSize + i]));
	}

	return logDet;
}