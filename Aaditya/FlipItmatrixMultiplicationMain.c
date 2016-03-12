// This is the main function of Matrix Multiplication function which will fault injection methods
// such as FLIPIT_Init(rank, argc, argv, seed); and FLIPIT_Finalize(NULL); 
// Fault is injected in multiplication function which is define matrixMultiplcationFunction
// In this example I have used Flipt LLVM based fault injection

//# build the executable
//gcc -I$FLIPIT_PATH/include -o main.o -c main.c
//gcc -o test final.o main.o -L$FLIPIT_PATH/lib/ -lcorrupt
//./test

#include <stdio.h>
#include "FlipIt/corrupt/corrupt.h"

void take_data(int a[][10], int b[][10], int r1,int c1, int r2, int c2);
void multiplication(int a[][10],int b[][10],int mult[][10],int r1,int c1,int r2,int c2);
void display(int mult[][10], int r1, int c2);
int main(int argc, char** argv)
{
    int a[10][10], b[10][10], mult[10][10], r1, c1, r2, c2, i, j, k;
    int seed = 7;
    int rank = 0;
    FLIPIT_Init(rank, argc, argv, seed);
    printf("Enter rows and column for first matrix: ");
    scanf("%d%d", &r1, &c1);
    printf("Enter rows and column for second matrix: ");
    scanf("%d%d",&r2, &c2);

/* If colum of first matrix in not equal to row of second matrix, asking user to enter the size of matrix again. */

    while (c1!=r2)
    {
        printf("Error! column of first matrix not equal to row of second.\n");
        printf("Enter rows and column for first matrix: ");
        scanf("%d%d", &r1, &c1);
        printf("Enter rows and column for second matrix: ");
        scanf("%d%d",&r2, &c2);
    }
    take_data(a,b,r1,c1,r2,c2);  /* Function to take matices data */
    multiplication(a,b,mult,r1,c1,r2,c2); /* Function to multiply two matrices. */
    display(mult,r1,c2); /* Function to display resultant matrix after multiplication. */
	FLIPIT_Finalize(NULL);    
	return 0;
}

void take_data(int a[][10], int b[][10], int r1,int c1, int r2, int c2)
{
    int i,j;
    printf("\nEnter elements of matrix 1:\n");
    for(i=0; i<r1; ++i)
    for(j=0; j<c1; ++j)
    {
        printf("Enter elements a%d%d: ",i+1,j+1);
        scanf("%d",&a[i][j]);
    }

    printf("\nEnter elements of matrix 2:\n");
    for(i=0; i<r2; ++i)
    for(j=0; j<c2; ++j)
    {
        printf("Enter elements b%d%d: ",i+1,j+1);
        scanf("%d",&b[i][j]);
    }
}



void display(int mult[][10], int r1, int c2)
{
    int i, j;
    printf("\nOutput Matrix:\n");
    for(i=0; i<r1; ++i)
    for(j=0; j<c2; ++j)
    {
        printf("%d  ",mult[i][j]);
        if(j==c2-1)
            printf("\n\n");
    }
}
