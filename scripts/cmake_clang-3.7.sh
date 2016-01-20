#!/usr/bin/env bash

path=$(cd >/dev/null 2>&1 $(dirname $0); pwd)
exec cmake \
	-GNinja \
	-DCUDA_HOST_COMPILER=$path/clang-legacy \
	-DCMAKE_CXX_COMPILER=/usr/local/clang+llvm-3.7.0-x86_64-apple-darwin/bin/clang++ \
	-DOpenMP_CXX_FLAGS='-fopenmp=libomp -I/usr/local/OpenMP-3.7.0-x86_64-apple-darwin14.4.0/include' \
	-DCMAKE_EXE_LINKER_FLAGS=-L/usr/local/OpenMP-3.7.0-x86_64-apple-darwin14.4.0/lib \
	"$@"
