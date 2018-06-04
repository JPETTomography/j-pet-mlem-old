/// \page test Testing
/// \brief More about unit-testing in this project
///
/// [catch]: https://github.com/philsquared/Catch
///
/// This project uses [Catch][catch] for unit-testing.
///
/// Tests are not build by default. In order to build and run tests run:
///
///     make test && ./test
///
/// Files containing tests have `_test.cpp` suffix.
///
/// Usage
/// -----
/// \verboutput test -h
///
/// \listoutput test -l

#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_USE_ANSI_COLOUR_CODES
#include "catch.hpp"
