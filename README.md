PET playground project
======================

This projects contains tools & test scripts for supporting PET project.

Authors
-------

* Piotr Bialas <pbialas@th.if.uj.edu.pl> (professor, project lead)
* Adam Strzelecki <adam.strzelecki@uj.edu.pl> (PhD student, code maintainer)
* Jakub Kowal <jakub.kowal@uj.edu.pl> (former PhD student and project member)

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
* [CUDA][cuda] 7.5 (automatically detected by `cmake`)
* [Ninja][ninja] 1.4 for faster re-builds (with `cmake -G Ninja`)
* `libpng` headers and libraries for PNG output
* [Boost][boost] 1.58 for advanced geometry calculation (`2d_barrel_geometry`)

[qtcreator]: http://qt-project.org/downloads
[cuda]: https://developer.nvidia.com/cuda-downloads
[ninja]: http://martine.github.io/ninja/
[libpng]: http://libpng.sourceforge.net
[boost]: http://www.boost.org

#### Notes

1. `libpng` is available on *Debian/Ubuntu GNU/Linux* via:

		apt-get install libpng12-dev

2. `libpng` is available on *Mac OS X* with [XQuartz][xquartz] X11 server

[xquartz]: http://xquartz.macosforge.org

3. `CUDA` 7.5 is required to compile GPU C++11 code

### Windows build requirements

[vs2015]: http://www.visualstudio.com/en-us/downloads/visual-studio-2015-downloads-vs

* [Visual Studio 2015][vs2015]
* *CMake* 3.0 for build script generation

##### Optional

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
using following methods:

1. Automatically using *Git* clean filter executed on each change introduced to
   index.

   This is set up running following command in root of the project:

       git config filter.format.clean $PWD/scripts/format-filter

2. Manually, calling following script in root of the project:

       ./scripts/format

3. Using *Qt Creator* 3.1 or newer *ClangFormat* beautifier.

   *File* format must be selected in *Settings > Beautifier > Clang Format >
   Style*.

When using *Qt Creator* code style used in this project can be imported using
*Settings > C++ > Code Style > Import* from `src/Google.xml` file.

Naming Convention
-----------------

1. Camel-case naming for classes and template names i.e. `SmallPotato`
2. Lower case with underscores for class files, i.e.: `small_potato.h`
3. Constants are upper case with underscores, i.e.: `BROWN_POTATO`
4. Variables and instances using lower case, i.e.: `some_potato`
5. Type template parameters should have `Type` suffix, i.e.: `SizeType`
6. Class template parameters should have `Class` suffix, i.e.: `DetectorClass`

Type Usage
----------

This project can be built to use `float` or `double` as floating point
representation, depending on the need, therefore:

1. All classes must be templates taking `FType` (aka floating point type),
   `SType` (aka short integer type) when necessary. They should not use `float`,
   `double` or `short` directly, eg.:

       template <typename FType> struct Point {
         using F = FType;
         F x, y
       };

2. All template instantiations shall use `F` and `S` types respectively
   including `common/types.h` prior instantiation, eg.:

       #include "2d/geometry/point.h"
       #include "2d/geometry/pixel.h"
       #include "common/types.h"
       
       using Point = PET2D::Point<F>;
       using Pixel = PET2D::Point<S>;

Serialization & Deserialization
-------------------------------

Project uses `std::istream` & `std::ostream` for textual
serialization/deserialization, and special `util::ibstream` & `util::obstream`
for binary serialization/deserialization.

Since this project is not using virtual methods, if class uses both textual and
binary serialization, it should define separate methods/constructors for both
`std::` text streams and `util::` binary streams.

#### Serialization

All classes that can be serialized into `ostream` or `obstream` should define
"friend" `operator<<` inside their body, eg.:

    std::ostream& operator<<(std::ostream& out, const Point& p) {
      return out << x << ' ' << y;
    }
    util::obstream& operator<<(util::obstream& out, const Point& p) {
      return out << x << y;  // NOTE: no ' ' separator since it is binary stream
    }

#### Deserialization

All classes that can be deserialized from input stream should provide
constructors taking input stream reference as an argument, eg.:

    Point(std::istream& in) { in >> x >> y; }
    Point(util::ibstream& in) { in >> x >> y; }

GPU (CUDA) Compatible Code
--------------------------

This project uses C++11 extensively, it also targets *NVIDIA CUDA* version 7.5
providing complete C++11 support.

Most of the classes as reused between CPU & GPU code using following rules:

1. All CUDA enabled class source files should include `util/cuda/compat.h`

2. All class methods available to CUDA should be started with single `_`, that
   will be turned into `__device__` when compiling via CUDA, eg.:

       _ Vector operator+(const Vector& other) { ... }

3. All CUDA enabled classes should NOT use `<cmath>` functions like `std::sin`,
   but instead use `compat::sin` counterparts which will be mapped to proper
   CUDA functions.

4. All CUDA enabled classes should NOT use any STL containers or headers, if it
   is necessary class may provide CPU only methods using preprocessor
   conditional blocks, eg.:

       #if !__CUDACC__
       #include <ostream>
       #endif
       
       ...
       
       #if !__CUDACC__
       std::ostream& operator<<(std::ostream& out, const Point& p) {
         return out << x << ' ' << y;
       }
       #endif
