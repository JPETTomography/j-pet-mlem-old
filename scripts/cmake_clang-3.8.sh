#!/usr/bin/env bash

path=$(cd >/dev/null 2>&1 $(dirname $0); pwd)
exec cmake \
	-GNinja \
	-DCUDA_HOST_COMPILER=$path/clang-legacy \
	-DCMAKE_CXX_COMPILER=/usr/local/clang+llvm-3.8.0-x86_64-apple-darwin/bin/clang++ \
	-DOpenMP_CXX_FLAGS=-fopenmp=libomp \
	-DCMAKE_EXE_LINKER_FLAGS="{-L,-Wl\,-rpath\,}/usr/local/clang+llvm-3.8.0-x86_64-apple-darwin/lib" \
	"$@"
