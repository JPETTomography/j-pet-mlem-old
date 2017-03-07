#!/usr/bin/env bash

path=$(cd >/dev/null 2>&1 $(dirname $0); pwd)
base=
for d in /usr/local/Cellar/llvm/* ; do
	[ -d $d ] && base="$d"
done

if [ ! -d $base ]; then
	echo >&2 'Cannot find Clang in /usr/local/Cellar/llvm'
	exit 1
fi

echo 'Using Clang at:' $base

exec cmake \
	-GNinja \
	-DPNG_PNG_INCLUDE_DIR=/usr/local/include \
	-DCUDA_HOST_COMPILER=$path/clang-legacy \
	-DCMAKE_CXX_COMPILER=$base/bin/clang++ \
	-DOpenMP_CXX_FLAGS=-fopenmp=libomp \
	-DCMAKE_EXE_LINKER_FLAGS="{-L,-Wl\,-rpath\,}$base/lib" \
	"$@"
