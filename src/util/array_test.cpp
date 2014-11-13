#include "test.h"
#include "array.h"

struct Foo {
  int foo;
  Foo(int foo) : foo(foo) {}
};

struct Bar {
  int foo, bar;
  Bar(int foo, int bar) : foo(foo), bar(bar) {}
};

struct Boo {
  int* boo;
  Boo() : boo(new int) {}
  ~Boo() { delete boo; }
};

TEST("util/array") {

  SECTION("default_alignment") {
    using Foos = util::array<3, Foo>;
    CHECK(Foos::value_size == sizeof(int));
    CHECK(Foos::storage_size == sizeof(int));
    CHECK(Foos::alignment == alignof(int));

    using Bars = util::array<3, Bar>;
    CHECK(Bars::value_size == 2 * sizeof(int));
    CHECK(Bars::storage_size == 2 * sizeof(int));
    CHECK(Bars::alignment == alignof(int));
  }

  SECTION("no_non_trivial_desctructors") {
#if SHOULD_FAIL
    // this will fail with: "no member named"
    using Boos = util::array<3, Boo>;
#endif
    // we must explicitly specify storage type and alignment
    using Boos = util::array<3, Boo, alignof(Boo), Boo>;
    CHECK(Boos::value_size == sizeof(Boo));
    CHECK(Boos::storage_size == sizeof(Boo));
    CHECK(Boos::alignment == alignof(Boo));
  }

  SECTION("custom_alignment") {
    using Foos = util::array<3, Foo, 2 * alignof(int)>;
    CHECK(Foos::value_size == sizeof(int));
    CHECK(Foos::storage_size == 2 * sizeof(int));
    CHECK(Foos::alignment == 2 * alignof(int));

    using Bars = util::array<3, Bar, 2 * alignof(int)>;
    CHECK(Bars::value_size == 2 * sizeof(int));
    CHECK(Bars::storage_size == 2 * sizeof(int));
    CHECK(Bars::alignment == 2 * alignof(int));
  }

  SECTION("push_back") {
    util::array<3, Foo> foos;
    CHECK(foos.size() == 0);

    foos.push_back(Foo(1));
    CHECK(foos.size() == 1);
    CHECK(foos[0].foo == 1);

    foos.push_back(Foo(2));
    CHECK(foos.size() == 2);
    CHECK(foos[1].foo == 2);
  }

  SECTION("emplace_back") {
    util::array<3, Foo> foos;
    CHECK(foos.size() == 0);

    foos.emplace_back(1);
    CHECK(foos.size() == 1);
    CHECK(foos[0].foo == 1);

    foos.emplace_back(2);
    CHECK(foos.size() == 2);
    CHECK(foos[1].foo == 2);
  }

  SECTION("in_place_constructor") {
    util::array<4, Foo> foos{ 1, 2, 3, 4 };
    CHECK(foos.size() == 4);
    CHECK(foos[0].foo == 1);
    CHECK(foos[1].foo == 2);
    CHECK(foos[2].foo == 3);
    CHECK(foos[3].foo == 4);
  }

  SECTION("in_place_constructor_custom_aligned") {
    util::array<4, Foo, 2 * alignof(int)> foos{ 1, 2, 3, 4 };
    CHECK(foos.size() == 4);
    CHECK(foos[0].foo == 1);
    CHECK(foos[1].foo == 2);
    CHECK(foos[2].foo == 3);
    CHECK(foos[3].foo == 4);
  }
}
