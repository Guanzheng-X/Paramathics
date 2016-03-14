CC=/Users/wuyue/LLFI/llvm/bin/clang++
LINKER=/Users/wuyue/LLFI/llvm/bin/llvm-link
OUTPUT=CG.ll
CFLAGS=-w -emit-llvm -fno-use-cxa-atexit -S
LINKER_FLAGS=-o $(OUTPUT) -S
SRCDIR_OBJS=main.ll math.ll sparse_matrix.ll 

build:$(SRCDIR_OBJS)
	$(LINKER) $(LINKER_FLAGS) $(SRCDIR_OBJS)

main.ll: main.cpp
	$(CC) $(CFLAGS) main.cpp

math.ll: math.cpp
	$(CC) $(CFLAGS) math.cpp

sparse_matrix.ll: sparse_matrix.cpp
	$(CC) $(CFLAGS) sparse_matrix.cpp

clean:
	rm -rf *.ll *.bc llfi* CG.ll

