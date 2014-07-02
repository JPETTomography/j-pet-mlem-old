PET playground project
======================

This projects contains tools & test scripts for supporting PET project.

Members
-------

* Piotr Bialas <pbialas@th.if.uj.edu.pl> (professor, project lead)
* Jakub Kowal <jakub.kowal@uj.edu.pl> (PhD student)
* Adam Strzelecki <adam.strzelecki@uj.edu.pl> (PhD student, code maintainer)

Prerequisites
-------------

### Minimum build requirements

* UNIX compatible build environment such as *Linux* or *Mac OS X*
* *C++11* compatible compiler i.e. *GCC* 4.6, *Clang* 3.2 or *ICC* 13
* *CMake* 2.8 with *GNU Make* 3.8 or *Ninja* 1.4
* `libpng` headers and libraries for PNG output

### Optional

* *QtCreator* 3.1 for project editing via `CMakeLists.txt`
* *CUDA* 6.0 (automatically detected by `cmake`)

Coding style
------------

[style]: http://dev.chromium.org/developers/coding-style
[clang-format]: http://clang.llvm.org/docs/ClangFormat.html 

This project follows *C++11* and [Chromium/Google source coding style][style].

This coding style is enforced using [clang-format][clang-format] reformat via:

	./scripts/format

When using *Qt Creator* code style used in this project can be imported using
*Settings > C++ > Code Style > Import* from `src/Google.xml` file.

Build
-----

Use *Qt Creator* to build source code or build from command line using `cmake`
with following commands:

	cmake . && make

Additional options:

1. To build outside of project directory use:

		cmake <path_to_project>

3. To use different than default compiler:

		cmake -DCMAKE_CXX_COMPILER=icpc  # for Intel C++ compiler
		cmake -DCMAKE_CXX_COMPILER=g++   # for GCC
		cmake -DCMAKE_CXX_COMPILER=clang # for Clang
