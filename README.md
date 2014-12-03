PET playground project
======================

This projects contains tools & test scripts for supporting PET project.

Authors
-------

* Piotr Bialas <pbialas@th.if.uj.edu.pl> (professor, project lead)
* Adam Strzelecki <adam.strzelecki@uj.edu.pl> (PhD student, code maintainer)
* Jakub Kowal <jakub.kowal@uj.edu.pl> (PhD student)

Prerequisites
-------------

### UNIX build requirements

* UNIX compatible build environment such as *Linux* or *Mac OS X*
* *C++11* compatible compiler i.e. *GCC* 4.6, *Clang* 3.2 or *ICC* 13
* [CMake][cmake] 2.8 for build script generation
* *GNU Make* 3.8 for building using `Makefile`

[cmake]: http://www.cmake.org

##### Optional

* [QtCreator][qtcreator] 3.1 for project editing via `CMakeLists.txt`
* [CUDA][cuda] 7.0 (automatically detected by `cmake`)
* [Ninja][ninja] 1.4 for faster re-builds (with `cmake -G Ninja`)
* `libpng` headers and libraries for PNG output

[qtcreator]: http://qt-project.org/downloads
[cuda]: https://developer.nvidia.com/cuda-downloads
[ninja]: http://martine.github.io/ninja/
[libpng]: http://libpng.sourceforge.net

#### Notes

1. `libpng` is available on *Debian/Ubuntu GNU/Linux* via:

		apt-get install libpng12-dev

2. `libpng` is available on *Mac OS X* with [XQuartz][xquartz] X11 server

[xquartz]: http://xquartz.macosforge.org

3. `CUDA` 7.0 is required to compile GPU C++11 code

### Windows build requirements

* *Visual Studio* 2013 (aka *MSVC* 12)
* *CMake* 3.0 for build script generation

##### Optional

* [CUDA][cuda] 6.5 (automatically detected by `cmake`)
* [CMake Tools][cmaketools] for editing *CMake* in *Visual Studio*

[cmaketools]: http://cmaketools.codeplex.com

#### Notes

1. Although project now compiles on *Windows*, however our main build
   environment is UNIX, so *Windows* compatibility may break any time.

2. *Windows* build flavor automatically downloads and compiles `libpng` and
   `zlib`, since these both does not come with *Windows* by default. Therefore
   build requires Internet connection.

3. *Ninja* currently does not work on *Windows* with our project.

Building
--------

Project uses *CMake* to generate platform specific build scripts, to build with
default settings run:

	cmake . && make

Additional options:

1. To build outside of project directory use:

		cmake <path_to_project>

2. To use *Ninja* build instead default (`make`) use:

		cmake -G Ninja

3. To use different than default compiler:

		cmake -DCMAKE_CXX_COMPILER=icpc  # for Intel C++ compiler
		cmake -DCMAKE_CXX_COMPILER=g++   # for GCC
		cmake -DCMAKE_CXX_COMPILER=clang # for Clang

Build Configuration
-------------------

There are several compile time options that can be changed, enabled or
disabled. All *PET* specific options are *CMake* variables/options prefixed
with `PET_`.

To list all available configure options:

	cmake -LH

Once command above is run, typing `cmake -DPET` and pressing **[TAB]** in
terminal completes all PET specific variables, unless *Bash* completion is
disabled.

To disable some option enabled by default, e.g. `PET_WARP_GRANULARITY`:

	cmake -DPET_WARP_GRANULARITY:BOOL=OFF

Testing
-------

Tests are not build by default. In order to build and run tests run:

	make test && ./test

Files containing tests have `_test.cpp` suffix.

Documentation
-------------

[doxygen]: http://www.doxygen.org/
[markdown]: http://daringfireball.net/projects/markdown/

Project uses [Doxygen][doxygen] to generate documentation out of source code
and [Markdown][markdown] files. To generate documentation run:

	doxygen

Next open `doc/html/index.html` file to view documentation.

Coding Style
------------

[style]: http://dev.chromium.org/developers/coding-style
[clang-format]: http://clang.llvm.org/docs/ClangFormat.html

This project follows *C++11* and [Chromium/Google source coding style][style]
with custom settings described in `.clang-format`.

Prior committing code should be formatted using [clang-format][clang-format]
script calling:

	./scripts/format

When using *Qt Creator* code style used in this project can be imported using
*Settings > C++ > Code Style > Import* from `src/Google.xml` file.

*Qt Creator* 3.1 supports also *ClangFormat* beautifier, however *File* format
should be selected in *Settings > Beautifier > Clang Format > Style*.

Naming Convention
-----------------

1. Camel-case naming for classes and template names i.e. `SmallPotato`
2. Lower case with underscores for class files, i.e.: `small_potato.h`
3. Constants are upper case with underscores, i.e.: `BROWN_POTATO`
4. Variables and instances using lower case, i.e.: `some_potato`
