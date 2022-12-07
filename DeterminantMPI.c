#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>


#define ARRAY_SIZE 16

double determinantOfMatrix(double* mat, int n);
double logDeterminantOfMatrix(double* mat, int n);


int main(int argc, char* argv[])
{
	char fileName[50];
	double* matrix = calloc(ARRAY_SIZE * ARRAY_SIZE, sizeof(double));

	int myRank, commSize;

	// Initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	if(myRank == 0)
	{
		// Create filename
		sprintf(fileName, "../DET_DATA/m0016x0016.bin");
		printf("Reading array file %s of size %dx%d\n", fileName, ARRAY_SIZE, ARRAY_SIZE);
		// Open file
		FILE* dataFile = fopen(fileName, "rb");
		// Read elements
		int i, j;
		for(i = 0; i < ARRAY_SIZE; i++)
		{
			for(j = 0; j < ARRAY_SIZE; j++)
			{
				fread(&matrix[i * ARRAY_SIZE + j], sizeof(double), 1, dataFile);
				printf("matrix[%d][%d]=%f\n", i, j, matrix[i * ARRAY_SIZE + j]);
			}
		}
		printf("Matrix has been read.\n");
	}

	//Broadcast the matrix to all processes
	MPI_Bcast(matrix, ARRAY_SIZE * ARRAY_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	printf("Rank: %d\t Example: %f\n", myRank, matrix[0]);

	int myStart = myRank * ARRAY_SIZE / commSize;
	int myEnd = (myRank + 1) * ARRAY_SIZE / commSize;

	double myDet = 0;
	double myLogDet = 0;
	// Loop for my determinants
	for(int i = myStart; i < myEnd; i++)
	{
		double currDet;
		double currLogDet;
		int newJ, newK;
		//copy the sub matrix
		double* subMatrix = calloc((ARRAY_SIZE - 1) * (ARRAY_SIZE - 1), sizeof(double));
		for(int j = 1; j < ARRAY_SIZE; j++)
		{
			for(int k = 0; k < ARRAY_SIZE; k++)
			{
				if(k == i)
				{
					continue;
				}
				newJ = j - 1;
				newK = k > i ? k - 1 : k;
				subMatrix[newJ * (ARRAY_SIZE - 1) + newK] = matrix[j * ARRAY_SIZE + k];
			}
		}

		// Calculate determinant
		currDet = matrix[i] * determinantOfMatrix(subMatrix, ARRAY_SIZE - 1);
		currLogDet = log10(fabs(matrix[i])) + logDeterminantOfMatrix(subMatrix, ARRAY_SIZE - 1);

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

	// Reduce and sum all determinants
	double globalDet;
	double globalLogDet;
	MPI_Reduce(&myDet, &globalDet, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&myLogDet, &globalLogDet, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	free(matrix);

	if(myRank == 0)
	{
		printf("Determinant of matrix is %f\n", globalDet);
		printf("Log Determinant of matrix is %f\n", globalLogDet);
	}

	MPI_Finalize();
}


double determinantOfMatrix(double* mat, int n)
{
	int i, j, k;
	double det = 1;
	double ratio;

	/* Applying Gauss Elimination */
	for(i = 0; i < n; i++)
	{
		if(mat[i * n + i] == 0.0)
		{
			printf("Mathematical Error!");
			exit(0);
		}
		for(j = i + 1; j < n; j++)
		{
			ratio = mat[j * n + i] / mat[i * n + i];

			for(k = 0; k < n; k++)
			{
				mat[j * n + k] = mat[j * n + k] - ratio * mat[i * n + k];
			}
		}
	}

	/* Finding determinant by multiplying
	 elements in principal diagonal elements */
	for(i = 0; i < n; i++)
	{
		det = det * mat[i * n + i];
	}

	return det;
}


double logDeterminantOfMatrix(double* mat, int n)
{
	int i, j, k;
	double det = 0;
	double ratio;

	/* Applying Gauss Elimination */
	for(i = 0; i < n; i++)
	{
		if(mat[i * n + i] == 0.0)
		{
			printf("Mathematical Error!");
			exit(0);
		}
		for(j = i + 1; j < n; j++)
		{
			ratio = mat[j * n + i] / mat[i * n + i];

			for(k = 0; k < n; k++)
			{
				mat[j * n + k] = mat[j * n + k] - ratio * mat[i * n + k];
			}
		}
	}

	/* Finding determinant by multiplying
	 elements in principal diagonal elements */
	for(i = 0; i < n; i++)
	{
		det = det + log10(fabs(mat[i * n + i]));
	}

	return det;
}
