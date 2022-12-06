#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define ARRAY_SIZE 16

int determinantOfMatrix(int* mat, int n);
void swap(int* p, int* q);


int main(int argc, char* argv[])
{
	char fileName[50];
	double a[ARRAY_SIZE][ARRAY_SIZE];
	int i, j;
	double det;
	// Create filename
	sprintf(fileName, "DET_DATA/m0032x0032.bin");
	printf("Reading array file %s of size %dx%d\n", fileName, ARRAY_SIZE, ARRAY_SIZE);
	// Open file
	FILE* dataFile = fopen(fileName, "rb");
	// Read elements
	for(i = 0; i < ARRAY_SIZE; i++)
	{
		for(j = 0; j < ARRAY_SIZE; j++)
		{
			fread(&a[i][j], sizeof(double), 1, dataFile);
			printf("a[%d][%d]=%f\n", i, j, a[i][j]);
		}
	}
	printf("Matrix has been read.\n");

	// Calculate determinant

}


int determinantOfMatrix(int* mat, int n)
{
	int num1, num2, det = 1, index,
			total = 1; // Initialize result

	// temporary array for storing row
	int temp[n + 1];

	// loop for traversing the diagonal elements
	for(int i = 0; i < n; i++)
	{
		index = i; // initialize the index

		// finding the index which has non zero value
		while(index < n && mat[index * ARRAY_SIZE + i] == 0)
		{
			index++;
		}
		if(index == n) // if there is non zero element
		{
			// the determinant of matrix as zero
			continue;
		}
		if(index != i)
		{
			// loop for swapping the diagonal element row and
			// index row
			for(int j = 0; j < n; j++)
			{
				swap(&mat[index * ARRAY_SIZE + j], &mat[i * ARRAY_SIZE + j]);
			}
			// determinant sign changes when we shift rows
			// go through determinant properties
			det = det * pow(-1, index - i);
		}

		// storing the values of diagonal row elements
		for(int j = 0; j < n; j++)
		{
			temp[j] = mat[i * ARRAY_SIZE + j];
		}
		// traversing every row below the diagonal element
		for(int j = i + 1; j < n; j++)
		{
			num1 = temp[i]; // value of diagonal element
			num2 = mat[j * ARRAY_SIZE + i]; // value of next row element

			// traversing every column of row
			// and multiplying to every row
			for(int k = 0; k < n; k++)
			{
				// multiplying to make the diagonal
				// element and next row element equal
				mat[j * ARRAY_SIZE + k]
						= (num1 * mat[j * ARRAY_SIZE + k]) - (num2 * temp[k]);
			}
			total = total * num1; // Det(kA)=kDet(A);
		}
	}

	// multiplying the diagonal elements to get determinant
	for(int i = 0; i < n; i++)
	{
		det = det * mat[i * ARRAY_SIZE + i];
	}
	return (det / total); // Det(kA)/k=Det(A);
}


void swap(int* p, int* q)
{
	//p=&n1 so p store the address of n1, so *p store the value of n1
	//q=&n2 so q store the address of n2, so *q store the value of n2

	int tmp;
	tmp = *p; // tmp store the value of n1
	*p = *q;    // *p store the value of *q that is value of n2
	*q = tmp;   // *q store the value of tmp that is the value of n1
}