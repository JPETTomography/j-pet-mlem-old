PET playground project
======================

This projects contains tools & test scripts for supporting PET project.

Members
-------

* Piotr Bialas <pbialas@th.if.uj.edu.pl> (head)
* Jakub Kowal <jakub.kowal@uj.edu.pl> (PhD student)
* Adam Strzelecki <adam.strzelecki@uj.edu.pl> (PhD student)

Prerequisites
-------------

### Build

* UNIX compatible environment (i.e. *Linux* or *Mac OS X*)
* One of *GCC 4.5* (or higher), *Clang 3.2* (or higher)

### Source code editing

* *QtCreator 4.6* (or higher) for convenient C++ editing source files via `stc/PET.pro` project file
* *TextMate 2.x*, project contains `.tm_properties`
* *Emacs 23* (or higher), project contains `.dir-locals`

Coding style
------------

This project follows C++11 and [Chromium/Google source coding
style](http://dev.chromium.org/developers/coding-style). This coding style is enforced using
[clang-format](http://clang.llvm.org/docs/ClangFormat.html) reformat via:

	make style

Build
-----

To build project file go into `src/` folder and type:

	make

To build with particular optimization (i.e. `0`) use:

	make O=0

To build with particular compiler (i.e. `gcc`) use:

	make CC=gcc

To build with *OpenMP* support use:

	make OMP=1

To display complete build commands use:

	make Q=
