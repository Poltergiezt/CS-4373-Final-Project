#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>


double determinantOfMatrix(double* mat, int matrixSize);
double logDeterminantOfMatrix(double* mat, int matrixSize);


int main(int argc, char* argv[])
{
	int myRank, commSize;
	// Initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	double myProgramStart = MPI_Wtime();

	double* matrix;
	int matrixSize;
	if(myRank == 0)
	{
		if(argc < 2)
		{
			fprintf(stderr, "Usage: %s matrix_size\n", argv[0]);
			MPI_Finalize();
			exit(1);
		}
		matrixSize = strtol(argv[1], NULL, 10);
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
			MPI_Finalize();
			exit(1);
		}
		printf("(1) %s\n", fileName);
		printf("(2) %d\n", matrixSize);
		// Read elements
		matrix = calloc(matrixSize * matrixSize, sizeof(double));
		for(int i = 0; i < matrixSize; i++)
		{
			for(int j = 0; j < matrixSize; j++)
			{
				fread(&matrix[i * matrixSize + j], sizeof(double), 1, dataFile);
			}
		}
		fclose(dataFile);
	}

	MPI_Bcast(&matrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(myRank != 0)
	{
		matrix = calloc(matrixSize * matrixSize, sizeof(double));
	}
	MPI_Bcast(matrix, matrixSize * matrixSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int myStart;
	int myEnd;
	if(matrixSize % commSize != 0)
	{
		if(myRank < matrixSize % commSize)
		{
			// Can implement the ceiling function by allowing truncation and then adding 1.
			myStart = myRank * (matrixSize / commSize + 1);
			myEnd = (myRank + 1) * (matrixSize / commSize + 1);
		}
		else
		{
			// Default behavior of integer division is to truncate, so no special logic needed here.
			myStart = myRank * (matrixSize / commSize) + matrixSize % commSize;
			myEnd = (myRank + 1) * (matrixSize / commSize) + matrixSize % commSize;
		}
	}
	else
	{
		myStart = myRank * (matrixSize / commSize);
		myEnd = (myRank + 1) * (matrixSize / commSize);
	}

	double myDet = 0;
	double myLogDet = 0;
	double myParallelStart = MPI_Wtime();
	// Loop for my determinants
	for(int i = myStart; i < myEnd; i++)
	{
		double currDet;
		double currLogDet;
		int newJ, newK;
		//copy the sub matrix
		double* subMatrix = calloc((matrixSize - 1) * (matrixSize - 1), sizeof(double));
		for(int j = 1; j < matrixSize; j++)
		{
			for(int k = 0; k < matrixSize; k++)
			{
				if(k == i)
				{
					continue;
				}
				newJ = j - 1;
				newK = k > i ? k - 1 : k;
				subMatrix[newJ * (matrixSize - 1) + newK] = matrix[j * matrixSize + k];
			}
		}

		// Calculate determinant
		currDet = matrix[i] * determinantOfMatrix(subMatrix, matrixSize - 1);
		currLogDet = log10(fabs(matrix[i])) + logDeterminantOfMatrix(subMatrix, matrixSize - 1);

		//Determine sign
		if(i % 2 != 0)
		{
			currDet = -currDet;
			currLogDet = -currLogDet;
		}

		// Add to total
		myDet += currDet;
		myLogDet += currLogDet;
		free(subMatrix);
	}

	free(matrix);

	// Reduce and sum all determinants.
	double det;
	double logDet;
	MPI_Reduce(&myDet, &det, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&myLogDet, &logDet, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// Collect timings.
	double myParallelEnd = MPI_Wtime();
	double myProgramEnd = MPI_Wtime();
	double myParallelTime = myParallelEnd - myParallelStart;
	double myProgramTime = myProgramEnd - myProgramStart;
	double parallelTime;
	MPI_Reduce(&myParallelTime, &parallelTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	double programTime;
	MPI_Reduce(&myProgramTime, &programTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if(myRank == 0)
	{
		printf("(3) %f\n", det);
		printf("(4) %f\n", logDet);
		printf("%f\n", parallelTime);
		printf("%f\n", programTime);
	}

	MPI_Finalize();
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
