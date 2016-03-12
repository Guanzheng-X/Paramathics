// Flipt will inject fault in this method
// Below is the example of how to compile it with corrupt library
//$LLVM_BUILD_PATH/bin/clang  -g -emit-llvm work.c  -c -o work.bc 
//$LLVM_BUILD_PATH/bin/llvm-link $FLIPIT_PATH/include/FlipIt/corrupt/corrupt.bc work.bc  -o crpt_work.bc
//$LLVM_BUILD_PATH/bin/opt -load $FLIPIT_PATH/lib/libFlipItPass.so -FlipIt -srcFile "work.c" -singleInj 1 -prob 0.95 -byte -1 -arith 1 -ctrl 0 -ptr 0 -funcList "" crpt_work.bc  -o final.bc
//$LLVM_BUILD_PATH/bin/clang  -g -c final.bc -o final.o 

void multiplication(int a[][10],int b[][10],int mult[][10],int r1,int c1,int r2,int c2)
{
    int i,j,k;
/* Initializing elements of matrix mult to 0.*/
    for(i=0; i<r1; ++i)
    for(j=0; j<c2; ++j)
    {
       mult[i][j]=0;
    }
/* Multiplying matrix a and b and storing in array mult. */
    for(i=0; i<r1; ++i)
    for(j=0; j<c2; ++j)
    for(k=0; k<c1; ++k)
    {
        mult[i][j]+=a[i][k]*b[k][j];
    }
}
