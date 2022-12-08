#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define ARRAY_SIZE 16

double determinantOfMatrix(double* mat, int n);
double logDeterminantOfMatrix(double* mat, int n);


int main(int argc, char* argv[])
{
	char fileName[50];
	double* matrix = calloc(ARRAY_SIZE * ARRAY_SIZE, sizeof(double));
	int i, j;
	double det;
	// Create filename
	sprintf(fileName, "../DET_DATA/m0016x0016.bin");
	printf("Reading array file %s of size %dx%d\n", fileName, ARRAY_SIZE, ARRAY_SIZE);
	// Open file
	FILE* datafile = fopen(fileName, "rb");
	// Read elements
	for(i = 0; i < ARRAY_SIZE; i++)
	{
		for(j = 0; j < ARRAY_SIZE; j++)
		{
			fread(&matrix[i * ARRAY_SIZE + j], sizeof(double), 1, datafile);
			printf("matrix[%d][%d]=%f\n", i, j, matrix[i * ARRAY_SIZE + j]);
		}
	}
	printf("Matrix has been read.\n");

	//Calculate determinant
	det = determinantOfMatrix(matrix, ARRAY_SIZE);
	printf("Determinant of the matrix is: %f\n", det);
	double logDet = logDeterminantOfMatrix(matrix, ARRAY_SIZE);
	printf("Log determinant of the matrix is: %f\n", logDet);
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