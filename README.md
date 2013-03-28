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

### Minimum build requirements

* UNIX compatible environment such as *Linux* or *Mac OS X*
* *GCC* 4.6, *Clang* 3.2 or *ICC* 13
* *GNU Make* 3.8
* `libpng` headers and libraries for PNG output

### Source code editing (one of these)

* *QtCreator 4.6* with C++ auto-completion via `stc/PET.pro`
* *TextMate 2.x* using `.tm_properties`
* *Emacs 23* using `.dir-locals`

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
