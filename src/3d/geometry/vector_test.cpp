#include "util/test.h"

#include "3d/geometry/vector.h"

using Vector = PET3D::Vector<float>;

TEST("3d/geometry/vector/init") {

  Vector vec(1.0f, 2.0f, 3.0f);

  CHECK(vec.x == 1.0_e7);
  CHECK(vec.y == 2.0_e7);
  CHECK(vec.z == 3.0_e7);
}

TEST("3d/geometry/vector/arithmetic_assignement") {
  {
    Vector vec1(1.0f, 2.0f, 3.0f);
    Vector vec2(0.10f, 0.2f, 0.3f);
    vec1 += vec2;

    CHECK(vec1.x == 1.1_e7);
    CHECK(vec1.y == 2.2_e7);
    CHECK(vec1.z == 3.3_e7);
  }
  {
    Vector vec1(1.0f, 2.0f, 3.0f);
    Vector vec2(0.10f, 0.2f, 0.3f);
    vec1 -= vec2;

    CHECK(vec1.x == 0.9_e7);
    CHECK(vec1.y == 1.8_e7);
    CHECK(vec1.z == 2.7_e7);
  }
  {
    Vector vec1(1.0f, 2.0f, 3.0f);
    vec1 *= 2.0f;

    CHECK(vec1.x == 2.0_e7);
    CHECK(vec1.y == 4.0_e7);
    CHECK(vec1.z == 6.0_e7);
  }
  {
    Vector vec1(1.0f, 2.0f, 3.0f);
    vec1 /= 2.0f;

    CHECK(vec1.x == 0.5_e7);
    CHECK(vec1.y == 1.0_e7);
    CHECK(vec1.z == 1.5_e7);
  }
}

TEST("3d/geometry/vector/arithmetics") {
  {
    Vector lhs(1.0f, 2.0f, 3.0f);
    Vector rhs(0.10f, 0.2f, 0.3f);
    Vector vec = lhs + rhs;

    CHECK(vec.x == 1.1_e7);
    CHECK(vec.y == 2.2_e7);
    CHECK(vec.z == 3.3_e7);
  }
  {
    Vector lhs(1.0f, 2.0f, 3.0f);
    Vector rhs(0.10f, 0.2f, 0.3f);
    Vector vec = lhs - rhs;

    CHECK(vec.x == 0.9_e7);
    CHECK(vec.y == 1.8_e7);
    CHECK(vec.z == 2.7_e7);
  }
}

TEST("3d/geometry/vector/logical") {
  {
    Vector lhs(1.0f, 2.0f, 3.0f);
    Vector rhs(1.0f, 2.0f, 3.0f);

    CHECK(lhs == rhs);
    CHECK(!(lhs != rhs));
  }
  {
    Vector lhs(1.0f, 2.0f, 4.0f);
    Vector rhs(1.0f, 2.0f, 3.0f);

    CHECK(!(lhs == rhs));
    CHECK((lhs != rhs));
  }
  {
    Vector lhs(1.0f, 5.0f, 3.0f);
    Vector rhs(1.0f, 2.0f, 3.0f);

    CHECK(!(lhs == rhs));
    CHECK((lhs != rhs));
  }
  {
    Vector lhs(5.0f, 2.0f, 3.0f);
    Vector rhs(1.0f, 2.0f, 3.0f);

    CHECK(!(lhs == rhs));
    CHECK((lhs != rhs));
  }
}
