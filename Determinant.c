#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define ARRAYSIZE 32

double determinantOfMatrix(double * mat, int n);
void swap(double *p,double *q);


int main(int argc, char *argv[])
{
    char f_name[50];
    double * a = calloc(ARRAYSIZE*ARRAYSIZE, sizeof(double));
    int i,j;
    double det;
    //Create filename
    sprintf(f_name,"../DET_DATA/m0032x0032.bin");
    printf("Reading array file %s of size %dx%d\n",f_name,ARRAYSIZE,ARRAYSIZE);
    //Open file
    FILE *datafile=fopen(f_name,"rb");
    //Read elelements
    for (i=0; i< ARRAYSIZE; i++)
        for (j=0; j< ARRAYSIZE; j++)
        {
            fread(&a[i * ARRAYSIZE +j],sizeof(double),1,datafile);
            printf("a[%d][%d]=%f\n",i,j,a[i * ARRAYSIZE +j]);
        }
    printf("Matrix has been read.\n");

    //Calculate determinant
    det = determinantOfMatrix(a,ARRAYSIZE);
    printf("Determinant of the matrix is: %f\n",det);
}

double determinantOfMatrix(double *mat, int n)
{
    double num1, num2, det = 1,
            total = 1; // Initialize result
    int  index;

    // temporary array for storing row
    double temp[n + 1];

    // loop for traversing the diagonal elements
    for (int i = 0; i < n; i++)
    {
        index = i; // initialize the index

        // finding the index which has non zero value
        while (index < n && mat[index * ARRAYSIZE + i] == 0)
        {
            index++;
        }
        if (index == n) // if there is non zero element
        {
            // the determinant of matrix as zero
            continue;
        }
        if (index != i)
        {
            // loop for swapping the diagonal element row and
            // index row
            for (int j = 0; j < n; j++)
            {
                swap(&mat[index *ARRAYSIZE + j], &mat[i * ARRAYSIZE + j]);
            }
            // determinant sign changes when we shift rows
            // go through determinant properties
            det = det * pow(-1, index - i);
        }

        // storing the values of diagonal row elements
        for (int j = 0; j < n; j++)
        {
            temp[j] = mat[i * ARRAYSIZE + j];
        }
        // traversing every row below the diagonal element
        for (int j = i + 1; j < n; j++)
        {
            num1 = temp[i]; // value of diagonal element
            num2 = mat[j * ARRAYSIZE + i]; // value of next row element

            // traversing every column of row
            // and multiplying to every row
            for (int k = 0; k < n; k++)
            {
                // multiplying to make the diagonal
                // element and next row element equal
                mat[j * ARRAYSIZE + k]
                        = (num1 * mat[j * ARRAYSIZE + k]) - (num2 * temp[k]);
            }
            total = total * num1; // Det(kA)=kDet(A);
        }
    }

    // multiplying the diagonal elements to get determinant
    for (int i = 0; i < n; i++)
    {
        det = det * mat[i * ARRAYSIZE + i];
    }
    return (det / total); // Det(kA)/k=Det(A);
}

void swap(double *p,double *q)
{
    //p=&n1 so p store the address of n1, so *p store the value of n1
    //q=&n2 so q store the address of n2, so *q store the value of n2

    double tmp;
    tmp = *p; // tmp store the value of n1
    *p=*q;    // *p store the value of *q that is value of n2
    *q=tmp;   // *q store the value of tmp that is the value of n1
}