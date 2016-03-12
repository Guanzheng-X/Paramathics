// Flipt will inject fault in this method
// Below is the example of how to compile it with corrupt library
//$LLVM_BUILD_PATH/bin/clang  -g -emit-llvm work.c  -c -o work.bc 
//$LLVM_BUILD_PATH/bin/llvm-link $FLIPIT_PATH/include/FlipIt/corrupt/corrupt.bc work.bc  -o crpt_work.bc
//$LLVM_BUILD_PATH/bin/opt -load $FLIPIT_PATH/lib/libFlipItPass.so -FlipIt -srcFile "work.c" -singleInj 1 -prob 0.95 -byte -1 -arith 1 -ctrl 0 -ptr 0 -funcList "" crpt_work.bc  -o final.bc
//$LLVM_BUILD_PATH/bin/clang  -g -c final.bc -o final.o 


int work(int a, int b)
{
	return a + b;
}
