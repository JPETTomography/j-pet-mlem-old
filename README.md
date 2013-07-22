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
* C++11 compatible compiler i.e. *GCC* 4.6, *Clang* 3.2 or *ICC* 13
* *Qmake* from *Qt* 4.x or 5.x SDK
* *GNU Make* 3.8
* `libpng` headers and libraries for PNG output

### Source code editing (one of these)

* *QtCreator 4.7* with C++ auto-completion via `PET.pro`
* *TextMate 2.x* using `.tm_properties`
* *Emacs 23* using `.dir-locals`

Coding style
------------

This project follows C++11 and [Chromium/Google source coding
style](http://dev.chromium.org/developers/coding-style). This coding style is enforced using
[clang-format](http://clang.llvm.org/docs/ClangFormat.html) reformat via:

	./scripts/format

Build
-----

Use *Qt Creator* to build source code or build from command line using `qmake` with following commands:

	qmake
	make
