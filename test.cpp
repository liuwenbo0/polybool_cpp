#include "Polybool.hpp"
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace polybool;

// ==========================================
// 1. Test helpers: tolerant comparisons and formatting
// ==========================================

const double TEST_EPSILON = 1e-6;

bool isFloatEq(double a, double b) { return std::abs(a - b) < TEST_EPSILON; }

bool isVec2Eq(const Vec2 &a, const Vec2 &b) {
  return isFloatEq(a.x, b.x) && isFloatEq(a.y, b.y);
}

// Format PolyPoint values for debugging.
std::string pointToStr(const PolyPoint &p) {
  std::ostringstream oss;
  if (std::holds_alternative<Vec2>(p)) {
    auto v = std::get<Vec2>(p);
    oss << "[" << v.x << ", " << v.y << "]";
  } else {
    auto v = std::get<Vec6>(p);
    oss << "[" << v.v[0] << ", " << v.v[1] << ", " << v.v[2] << ", " << v.v[3]
        << ", " << v.v[4] << ", " << v.v[5] << "]";
  }
  return oss.str();
}

std::string regionToStr(const std::vector<PolyPoint> &region) {
  std::ostringstream oss;
  oss << "{";
  for (size_t i = 0; i < region.size(); ++i) {
    if (i > 0)
      oss << ", ";
    oss << pointToStr(region[i]);
  }
  oss << "}";
  return oss.str();
}

// Assert that two polygons are exactly equal.
void assertEqual(const Polygon &a, const Polygon &b) {
  if (a.inverted != b.inverted)
    throw std::runtime_error("Inverted flags do not match");
  if (a.regions.size() != b.regions.size())
    throw std::runtime_error("Regions count does not match");

  for (size_t i = 0; i < a.regions.size(); ++i) {
    if (a.regions[i].size() != b.regions[i].size()) {
      throw std::runtime_error("Region point count does not match at index " +
                               std::to_string(i));
    }
    for (size_t j = 0; j < a.regions[i].size(); ++j) {
      const auto &pa = a.regions[i][j];
      const auto &pb = b.regions[i][j];

      if (pa.index() != pb.index())
        throw std::runtime_error("Point types do not match");

      if (std::holds_alternative<Vec2>(pa)) {
        auto va = std::get<Vec2>(pa), vb = std::get<Vec2>(pb);
        if (!isFloatEq(va.x, vb.x) || !isFloatEq(va.y, vb.y)) {
          throw std::runtime_error("Vec2 mismatch: " + pointToStr(pa) +
                                   " != " + pointToStr(pb));
        }
      } else {
        auto va = std::get<Vec6>(pa), vb = std::get<Vec6>(pb);
        for (int k = 0; k < 6; ++k) {
          if (!isFloatEq(va.v[k], vb.v[k])) {
            throw std::runtime_error("Vec6 mismatch: " + pointToStr(pa) +
                                     " != " + pointToStr(pb));
          }
        }
      }
    }
  }
}

std::vector<Vec2> toNormalizedVec2Region(const std::vector<PolyPoint> &region) {
  std::vector<Vec2> points;
  for (const auto &point : region) {
    if (!std::holds_alternative<Vec2>(point)) {
      throw std::runtime_error(
          "Expected a line-only polygon region, got curve point " +
          pointToStr(point));
    }

    Vec2 vec = std::get<Vec2>(point);
    if (points.empty() || !isVec2Eq(points.back(), vec)) {
      points.push_back(vec);
    }
  }

  if (points.size() > 1 && isVec2Eq(points.front(), points.back())) {
    points.pop_back();
  }

  bool changed = true;
  while (changed && points.size() >= 3) {
    changed = false;
    for (size_t i = 0; i < points.size(); ++i) {
      const Vec2 &prev = points[(i + points.size() - 1) % points.size()];
      const Vec2 &curr = points[i];
      const Vec2 &next = points[(i + 1) % points.size()];
      double cross = (curr.x - prev.x) * (next.y - curr.y) -
                     (curr.y - prev.y) * (next.x - curr.x);
      if (std::abs(cross) < TEST_EPSILON) {
        points.erase(points.begin() + static_cast<long>(i));
        changed = true;
        break;
      }
    }
  }

  return points;
}

bool isCyclicVec2RegionEq(const std::vector<PolyPoint> &aRegion,
                          const std::vector<PolyPoint> &bRegion) {
  std::vector<Vec2> a = toNormalizedVec2Region(aRegion);
  std::vector<Vec2> b = toNormalizedVec2Region(bRegion);
  if (a.size() != b.size())
    return false;
  if (a.empty())
    return true;

  for (size_t offset = 0; offset < b.size(); ++offset) {
    bool sameDirection = true;
    bool reversed = true;
    for (size_t i = 0; i < a.size(); ++i) {
      if (!isVec2Eq(a[i], b[(offset + i) % b.size()])) {
        sameDirection = false;
      }
      if (!isVec2Eq(a[i], b[(offset + b.size() - i) % b.size()])) {
        reversed = false;
      }
    }
    if (sameDirection || reversed)
      return true;
  }

  return false;
}

// Assert that line-only polygons are geometrically equivalent:
// Region order, start point, winding direction, and redundant collinear points may differ.
void assertEquivalentVec2Polygon(const Polygon &actual,
                                 const Polygon &expected) {
  if (actual.inverted != expected.inverted)
    throw std::runtime_error("Inverted flags do not match");
  if (actual.regions.size() != expected.regions.size()) {
    throw std::runtime_error("Regions count does not match: " +
                             std::to_string(actual.regions.size()) + " vs " +
                             std::to_string(expected.regions.size()));
  }

  std::vector<bool> matched(actual.regions.size(), false);
  for (const auto &expectedRegion : expected.regions) {
    bool found = false;
    for (size_t i = 0; i < actual.regions.size(); ++i) {
      if (!matched[i] &&
          isCyclicVec2RegionEq(actual.regions[i], expectedRegion)) {
        matched[i] = true;
        found = true;
        break;
      }
    }
    if (!found) {
      throw std::runtime_error("Could not find expected region " +
                               regionToStr(expectedRegion));
    }
  }
}

Polygon rect(double minX, double minY, double maxX, double maxY) {
  return {{{Vec2{minX, minY}, Vec2{maxX, minY}, Vec2{maxX, maxY},
            Vec2{minX, maxY}}},
          false};
}

Polygon polygon(std::initializer_list<std::vector<PolyPoint>> regions,
                bool inverted = false) {
  return {std::vector<std::vector<PolyPoint>>(regions), inverted};
}

Polygon emptyPolygon(bool inverted = false) { return {{}, inverted}; }

Polygon invertedPolygon(Polygon poly) {
  poly.inverted = true;
  return poly;
}

// ==========================================
// 2. Receiver and log assertion helpers
// ==========================================

struct LogItem {
  std::string op;
  std::vector<double> args;
};

class Receiver : public IPolyBoolReceiver {
public:
  std::vector<LogItem> log;

  void beginPath() override { log.push_back({"beginPath", {}}); }
  void moveTo(double x, double y) override {
    log.push_back({"moveTo", {x, y}});
  }
  void lineTo(double x, double y) override {
    log.push_back({"lineTo", {x, y}});
  }
  void bezierCurveTo(double cp1x, double cp1y, double cp2x, double cp2y,
                     double x, double y) override {
    log.push_back({"bezierCurveTo", {cp1x, cp1y, cp2x, cp2y, x, y}});
  }
  void closePath() override { log.push_back({"closePath", {}}); }

  std::vector<LogItem> done() { return log; }
};

void assertLogEqual(const std::vector<LogItem> &a,
                    const std::vector<LogItem> &b) {
  if (a.size() != b.size())
    throw std::runtime_error("Log size mismatch: " + std::to_string(a.size()) +
                             " vs " + std::to_string(b.size()));
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i].op != b[i].op)
      throw std::runtime_error("Log OP mismatch at " + std::to_string(i));
    if (a[i].args.size() != b[i].args.size())
      throw std::runtime_error("Log args size mismatch at " +
                               std::to_string(i));
    for (size_t j = 0; j < a[i].args.size(); ++j) {
      if (!isFloatEq(a[i].args[j], b[i].args[j]))
        throw std::runtime_error("Log arg mismatch at " + std::to_string(i));
    }
  }
}

// ==========================================
// 3. Test data
// ==========================================

Polygon triangle1 = {{{Vec2{0, 0}, Vec2{5, 10}, Vec2{10, 0}}}, false};

Polygon triangle2 = {{{Vec2{5, 0}, Vec2{10, 10}, Vec2{15, 0}}}, false};

Polygon box1 = {{{Vec2{0, 0}, Vec2{5, 0}, Vec2{5, -5}, Vec2{0, -5}}}, false};

Polygon curve1 = {{{Vec2{0, 0}, Vec6{0, -5, 10, -5, 10, 0}}}, false};

Polygon squareA = rect(0, 0, 10, 10);

Polygon squareOverlapB = rect(5, 5, 15, 15);

Polygon squareDisjointB = rect(20, 20, 30, 30);

Polygon squareEdgeTouchB = rect(10, 0, 20, 10);

Polygon squarePointTouchB = rect(10, 10, 20, 20);

Polygon squareInsideB = rect(2, 2, 8, 8);

Polygon twoSquares = {{{Vec2{0, 0}, Vec2{10, 0}, Vec2{10, 10}, Vec2{0, 10}},
                       {Vec2{20, 0}, Vec2{30, 0}, Vec2{30, 10},
                        Vec2{20, 10}}},
                      false};

Polygon horizontalCutter = rect(5, 0, 25, 10);

// ==========================================
// 4. Run tests
// ==========================================

int main() {
  GeometryEpsilon geo;
  PolyBool pb(&geo);

  int pass = 0;
  int fail = 0;

  auto runTest = [&](const std::string &name,
                     const std::function<void()> &func) {
    try {
      func();
      std::cout << "pass   " << name << "\n";
      pass++;
    } catch (const std::exception &err) {
      std::cout << "FAIL   " << name << "\n";
      std::cerr << name << " failed: " << err.what() << "\n";
      fail++;
    }
  };

  // Test: two triangles overlap in a smaller triangle.
  // Input: triangle1 intersect triangle2.
  // Expected: one triangle with points (10,0), (5,0), (7.5,5).
  runTest("triangles / intersect", [&]() {
    Polygon expected =
        polygon({{Vec2{10, 0}, Vec2{5, 0}, Vec2{7.5, 5}}});
    assertEqual(pb.intersect(triangle1, triangle2), expected);
  });

  // Test: two triangles overlap and should form a single outer boundary.
  // Input: triangle1 union triangle2.
  // Expected: one five-point polygon around both triangles.
  runTest("triangles / union", [&]() {
    Polygon expected = polygon({{Vec2{10, 10}, Vec2{7.5, 5}, Vec2{5, 10},
                                 Vec2{0, 0}, Vec2{15, 0}}});
    assertEqual(pb.union_(triangle1, triangle2), expected);
  });

  // Test: a rectangle and a cubic curve can be unioned.
  // Input: box1 union curve1.
  // Expected: one region preserving the clipped cubic Bezier segment.
  runTest("curve / union with box", [&]() {
    Polygon expected =
        polygon({{Vec2{10, 0}, Vec6{10, -2.5, 7.5, -3.75, 5, -3.75},
                  Vec2{5, -5}, Vec2{0, -5}, Vec2{0, 0}}});
    assertEqual(pb.union_(box1, curve1), expected);
  });

  // Test: overlapping rectangles should intersect in their shared area.
  // Input: squareA intersect squareOverlapB.
  // Expected: rectangle (5,5)-(10,10).
  runTest("rect overlap / intersect", [&]() {
    Polygon expected = rect(5, 5, 10, 10);
    assertEquivalentVec2Polygon(pb.intersect(squareA, squareOverlapB),
                                expected);
  });

  // Test: overlapping rectangles should union into one L-shaped boundary.
  // Input: squareA union squareOverlapB.
  // Expected: one eight-point polygon around both rectangles.
  runTest("rect overlap / union", [&]() {
    Polygon expected =
        polygon({{Vec2{0, 0}, Vec2{10, 0}, Vec2{10, 5}, Vec2{15, 5},
                  Vec2{15, 15}, Vec2{5, 15}, Vec2{5, 10}, Vec2{0, 10}}});
    assertEquivalentVec2Polygon(pb.union_(squareA, squareOverlapB), expected);
  });

  // Test: subtracting the second rectangle removes the top-right overlap.
  // Input: squareA difference squareOverlapB.
  // Expected: one L-shaped polygon covering squareA minus (5,5)-(10,10).
  runTest("rect overlap / difference", [&]() {
    Polygon expected =
        polygon({{Vec2{0, 0}, Vec2{10, 0}, Vec2{10, 5}, Vec2{5, 5},
                  Vec2{5, 10}, Vec2{0, 10}}});
    assertEquivalentVec2Polygon(pb.difference(squareA, squareOverlapB),
                                expected);
  });

  // Test: reverse difference removes the bottom-left overlap from B.
  // Input: squareA differenceRev squareOverlapB.
  // Expected: one L-shaped polygon covering squareOverlapB minus squareA.
  runTest("rect overlap / differenceRev", [&]() {
    Polygon expected =
        polygon({{Vec2{5, 10}, Vec2{10, 10}, Vec2{10, 5}, Vec2{15, 5},
                  Vec2{15, 15}, Vec2{5, 15}}});
    assertEquivalentVec2Polygon(pb.differenceRev(squareA, squareOverlapB),
                                expected);
  });

  // Test: xor keeps only the non-overlapping parts of both rectangles.
  // Input: squareA xor squareOverlapB.
  // Expected: two L-shaped regions, one from each rectangle.
  runTest("rect overlap / xor", [&]() {
    Polygon expected =
        polygon({{Vec2{0, 0}, Vec2{10, 0}, Vec2{10, 5}, Vec2{5, 5},
                  Vec2{5, 10}, Vec2{0, 10}},
                 {Vec2{5, 10}, Vec2{10, 10}, Vec2{10, 5}, Vec2{15, 5},
                  Vec2{15, 15}, Vec2{5, 15}}});
    assertEquivalentVec2Polygon(pb.xor_(squareA, squareOverlapB), expected);
  });

  // Test: disjoint rectangles have no shared area.
  // Input: squareA intersect squareDisjointB.
  // Expected: empty polygon.
  runTest("disjoint rects / intersect", [&]() {
    assertEquivalentVec2Polygon(pb.intersect(squareA, squareDisjointB),
                                emptyPolygon());
  });

  // Test: disjoint rectangles should remain two separate regions after union.
  // Input: squareA union squareDisjointB.
  // Expected: squareA and squareDisjointB as two regions.
  runTest("disjoint rects / union", [&]() {
    Polygon expected =
        polygon({{Vec2{0, 0}, Vec2{10, 0}, Vec2{10, 10}, Vec2{0, 10}},
                 {Vec2{20, 20}, Vec2{30, 20}, Vec2{30, 30},
                  Vec2{20, 30}}});
    assertEquivalentVec2Polygon(pb.union_(squareA, squareDisjointB), expected);
  });

  // Test: subtracting a disjoint polygon should not change the first polygon.
  // Input: squareA difference squareDisjointB.
  // Expected: squareA unchanged.
  runTest("disjoint rects / difference", [&]() {
    assertEquivalentVec2Polygon(pb.difference(squareA, squareDisjointB),
                                squareA);
  });

  // Test: union of identical polygons should be the same polygon.
  // Input: squareA union squareA.
  // Expected: squareA unchanged.
  runTest("identical rects / union", [&]() {
    assertEquivalentVec2Polygon(pb.union_(squareA, squareA), squareA);
  });

  // Test: intersection of identical polygons should be the same polygon.
  // Input: squareA intersect squareA.
  // Expected: squareA unchanged.
  runTest("identical rects / intersect", [&]() {
    assertEquivalentVec2Polygon(pb.intersect(squareA, squareA), squareA);
  });

  // Test: subtracting a polygon from itself should remove everything.
  // Input: squareA difference squareA.
  // Expected: empty polygon.
  runTest("identical rects / difference", [&]() {
    assertEquivalentVec2Polygon(pb.difference(squareA, squareA),
                                emptyPolygon());
  });

  // Test: xor of identical polygons should remove everything.
  // Input: squareA xor squareA.
  // Expected: empty polygon.
  runTest("identical rects / xor", [&]() {
    assertEquivalentVec2Polygon(pb.xor_(squareA, squareA), emptyPolygon());
  });

  // Test: rectangles sharing a full edge should merge into one rectangle.
  // Input: squareA union squareEdgeTouchB.
  // Expected: rectangle (0,0)-(20,10).
  runTest("shared edge rects / union", [&]() {
    Polygon expected = rect(0, 0, 20, 10);
    assertEquivalentVec2Polygon(pb.union_(squareA, squareEdgeTouchB),
                                expected);
  });

  // Test: shared boundary edge has no area, so intersection is empty.
  // Input: squareA intersect squareEdgeTouchB.
  // Expected: empty polygon.
  runTest("shared edge rects / intersect", [&]() {
    assertEquivalentVec2Polygon(pb.intersect(squareA, squareEdgeTouchB),
                                emptyPolygon());
  });

  // Test: xor of edge-touching rectangles equals their union because overlap has
  // no area.
  // Input: squareA xor squareEdgeTouchB.
  // Expected: rectangle (0,0)-(20,10).
  runTest("shared edge rects / xor", [&]() {
    Polygon expected = rect(0, 0, 20, 10);
    assertEquivalentVec2Polygon(pb.xor_(squareA, squareEdgeTouchB), expected);
  });

  // Test: rectangles touching at one point should not merge into one region.
  // Input: squareA union squarePointTouchB.
  // Expected: two square regions sharing only point (10,10).
  runTest("point-touch rects / union", [&]() {
    Polygon expected =
        polygon({{Vec2{0, 0}, Vec2{10, 0}, Vec2{10, 10}, Vec2{0, 10}},
                 {Vec2{10, 10}, Vec2{20, 10}, Vec2{20, 20},
                  Vec2{10, 20}}});
    assertEquivalentVec2Polygon(pb.union_(squareA, squarePointTouchB),
                                expected);
  });

  // Test: a single shared point has no area, so intersection is empty.
  // Input: squareA intersect squarePointTouchB.
  // Expected: empty polygon.
  runTest("point-touch rects / intersect", [&]() {
    assertEquivalentVec2Polygon(pb.intersect(squareA, squarePointTouchB),
                                emptyPolygon());
  });

  // Test: union with a contained polygon should keep only the outer boundary.
  // Input: squareA union squareInsideB.
  // Expected: squareA unchanged.
  runTest("contained rect / union", [&]() {
    assertEquivalentVec2Polygon(pb.union_(squareA, squareInsideB), squareA);
  });

  // Test: intersection with a contained polygon should return the inner polygon.
  // Input: squareA intersect squareInsideB.
  // Expected: squareInsideB.
  runTest("contained rect / intersect", [&]() {
    assertEquivalentVec2Polygon(pb.intersect(squareA, squareInsideB),
                                squareInsideB);
  });

  // Test: subtracting an inner rectangle should leave an outer boundary and a
  // hole boundary.
  // Input: squareA difference squareInsideB.
  // Expected: two regions: outer squareA and inner squareInsideB boundary.
  runTest("contained rect / difference", [&]() {
    Polygon expected =
        polygon({{Vec2{0, 0}, Vec2{10, 0}, Vec2{10, 10}, Vec2{0, 10}},
                 {Vec2{2, 2}, Vec2{8, 2}, Vec2{8, 8}, Vec2{2, 8}}});
    assertEquivalentVec2Polygon(pb.difference(squareA, squareInsideB),
                                expected);
  });

  // Test: xor with a contained rectangle equals outer minus inner.
  // Input: squareA xor squareInsideB.
  // Expected: two regions: outer squareA and inner squareInsideB boundary.
  runTest("contained rect / xor", [&]() {
    Polygon expected =
        polygon({{Vec2{0, 0}, Vec2{10, 0}, Vec2{10, 10}, Vec2{0, 10}},
                 {Vec2{2, 2}, Vec2{8, 2}, Vec2{8, 8}, Vec2{2, 8}}});
    assertEquivalentVec2Polygon(pb.xor_(squareA, squareInsideB), expected);
  });

  // Test: reverse difference of a contained polygon removes the inner polygon.
  // Input: squareA differenceRev squareInsideB.
  // Expected: empty polygon.
  runTest("contained rect / differenceRev", [&]() {
    assertEquivalentVec2Polygon(pb.differenceRev(squareA, squareInsideB),
                                emptyPolygon());
  });

  // Test: one cutter can intersect multiple input regions independently.
  // Input: twoSquares intersect horizontalCutter.
  // Expected: right half of first square and left half of second square.
  runTest("multi-region / intersect", [&]() {
    Polygon expected =
        polygon({{Vec2{5, 0}, Vec2{10, 0}, Vec2{10, 10}, Vec2{5, 10}},
                 {Vec2{20, 0}, Vec2{25, 0}, Vec2{25, 10},
                  Vec2{20, 10}}});
    assertEquivalentVec2Polygon(pb.intersect(twoSquares, horizontalCutter),
                                expected);
  });

  // Test: one cutter can subtract from multiple input regions independently.
  // Input: twoSquares difference horizontalCutter.
  // Expected: left strip of first square and right strip of second square.
  runTest("multi-region / difference", [&]() {
    Polygon expected =
        polygon({{Vec2{0, 0}, Vec2{5, 0}, Vec2{5, 10}, Vec2{0, 10}},
                 {Vec2{25, 0}, Vec2{30, 0}, Vec2{30, 10},
                  Vec2{25, 10}}});
    assertEquivalentVec2Polygon(pb.difference(twoSquares, horizontalCutter),
                                expected);
  });

  // Test: union with an inverted first polygon keeps the complement flag.
  // Input: inverted(squareA) union squareOverlapB.
  // Expected: inverted L-shaped boundary equal to squareA minus overlap.
  runTest("inverted first polygon / union", [&]() {
    Polygon expected =
        polygon({{Vec2{0, 0}, Vec2{10, 0}, Vec2{10, 5}, Vec2{5, 5},
                  Vec2{5, 10}, Vec2{0, 10}}},
                true);
    assertEquivalentVec2Polygon(pb.union_(invertedPolygon(squareA),
                                          squareOverlapB),
                                expected);
  });

  // Test: intersection with an inverted first polygon equals B minus A.
  // Input: inverted(squareA) intersect squareOverlapB.
  // Expected: non-inverted L-shaped part of squareOverlapB outside squareA.
  runTest("inverted first polygon / intersect", [&]() {
    Polygon expected =
        polygon({{Vec2{5, 10}, Vec2{10, 10}, Vec2{10, 5}, Vec2{15, 5},
                  Vec2{15, 15}, Vec2{5, 15}}});
    assertEquivalentVec2Polygon(pb.intersect(invertedPolygon(squareA),
                                             squareOverlapB),
                                expected);
  });

  // Test: difference with an inverted first polygon remains inverted.
  // Input: inverted(squareA) difference squareOverlapB.
  // Expected: inverted union boundary of squareA and squareOverlapB.
  runTest("inverted first polygon / difference", [&]() {
    Polygon expected =
        polygon({{Vec2{0, 0}, Vec2{10, 0}, Vec2{10, 5}, Vec2{15, 5},
                  Vec2{15, 15}, Vec2{5, 15}, Vec2{5, 10}, Vec2{0, 10}}},
                true);
    assertEquivalentVec2Polygon(pb.difference(invertedPolygon(squareA),
                                              squareOverlapB),
                                expected);
  });

  // Test: xor with an inverted first polygon toggles the inverted flag and keeps
  // both non-overlap boundaries.
  // Input: inverted(squareA) xor squareOverlapB.
  // Expected: inverted polygon with two L-shaped regions.
  runTest("inverted first polygon / xor", [&]() {
    Polygon expected =
        polygon({{Vec2{0, 0}, Vec2{10, 0}, Vec2{10, 5}, Vec2{5, 5},
                  Vec2{5, 10}, Vec2{0, 10}},
                 {Vec2{5, 10}, Vec2{10, 10}, Vec2{10, 5}, Vec2{15, 5},
                  Vec2{15, 15}, Vec2{5, 15}}},
                true);
    assertEquivalentVec2Polygon(pb.xor_(invertedPolygon(squareA),
                                         squareOverlapB),
                                expected);
  });

  // Test: Shape path builder and combine().intersect().output() produce the
  // documented path command sequence.
  // Input: two multi-region Shape objects built with moveTo/lineTo/closePath.
  // Expected: three closed output regions with the exact listed path commands.
  runTest("shape output / intersect example", [&]() {
    Receiver recv;

    // C++ chaining and object lifetime:
    // C++ methods return object references, so local variables are the safest
    // way to hold temporary shapes that will be combined.
    Shape shape1 = pb.shape();
    shape1.beginPath()
        .moveTo(50, 50)
        .lineTo(150, 150)
        .lineTo(190, 50)
        .closePath()
        .moveTo(130, 50)
        .lineTo(290, 150)
        .lineTo(290, 50)
        .closePath();

    Shape shape2 = pb.shape();
    shape2.beginPath()
        .moveTo(110, 20)
        .lineTo(110, 110)
        .lineTo(20, 20)
        .closePath()
        .moveTo(130, 170)
        .lineTo(130, 20)
        .lineTo(260, 20)
        .lineTo(260, 170)
        .closePath();

    shape1.combine(shape2).intersect().output(recv);

    std::vector<LogItem> expectedLog = {
        {"beginPath", {}},         {"moveTo", {110, 110}},
        {"lineTo", {50, 50}},      {"lineTo", {110, 50}},
        {"lineTo", {110, 110}},    {"closePath", {}},
        {"moveTo", {150, 150}},    {"lineTo", {178, 80}},
        {"lineTo", {130, 50}},     {"lineTo", {130, 130}},
        {"lineTo", {150, 150}},    {"closePath", {}},
        {"moveTo", {260, 131.25}}, {"lineTo", {178, 80}},
        {"lineTo", {190, 50}},     {"lineTo", {260, 50}},
        {"lineTo", {260, 131.25}}, {"closePath", {}}};

    assertLogEqual(recv.done(), expectedLog);
  });

  // Test: Shape transforms are applied before path output.
  // Input: triangle path with setTransform(3,0,0,2,100,200).
  // Expected: output coordinates transformed by x'=3x+100, y'=2y+200.
  runTest("shape output / transform", [&]() {
    Receiver recv;
    Shape shape = pb.shape();

    shape.setTransform(3, 0, 0, 2, 100, 200)
        .beginPath()
        .moveTo(50, 50)
        .lineTo(-10, 50)
        .lineTo(10, 10)
        .closePath()
        .output(recv);

    std::vector<LogItem> expectedLog = {
        {"beginPath", {}},      {"moveTo", {250, 300}}, {"lineTo", {70, 300}},
        {"lineTo", {130, 220}}, {"lineTo", {250, 300}}, {"closePath", {}}};

    assertLogEqual(recv.done(), expectedLog);
  });

  std::cout << "\nPass: " << pass << "\nFail: " << fail << "\n";
  if (fail > 0) {
    return 1; // Return non-zero to indicate failure in standard C++ practices
  }

  return 0;
}
