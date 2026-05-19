#include "Polybool.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <optional>

namespace polybool {

// Generate a portable high-precision PI value.
const double PI = std::acos(-1.0);

Roots3 GeometryEpsilon::solveCubicNormalized(double a, double b,
                                             double c) const {
  Roots3 res;
  double a3 = a / 3.0;
  double b3 = b / 3.0;
  double Q = a3 * a3 - b3;
  double R = a3 * (a3 * a3 - b / 2.0) + c / 2.0;

  if (std::abs(R) < epsilon && std::abs(Q) < epsilon) {
    res.push_back(-a3);
    return res;
  }

  double F = a3 * (a3 * (4.0 * a3 * c - b3 * b) - 2.0 * b * c) +
             4.0 * b3 * b3 * b3 + c * c;

  if (std::abs(F) < epsilon) {
    double sqrtQ = std::sqrt(Q);
    if (R > 0) {
      res.push_back(-2.0 * sqrtQ - a3);
      res.push_back(sqrtQ - a3);
    } else {
      res.push_back(-sqrtQ - a3);
      res.push_back(2.0 * sqrtQ - a3);
    }
    return res;
  }

  double Q3 = Q * Q * Q;
  double R2 = R * R;

  if (R2 < Q3) {
    double ratio = (R < 0 ? -1.0 : 1.0) * std::sqrt(R2 / Q3);
    if (ratio < -1.0)
      ratio = -1.0;
    if (ratio > 1.0)
      ratio = 1.0;

    double theta = std::acos(ratio);
    double norm = -2.0 * std::sqrt(Q);
    res.push_back(norm * std::cos(theta / 3.0) - a3);
    res.push_back(norm * std::cos((theta + 2.0 * PI) / 3.0) - a3);
    res.push_back(norm * std::cos((theta - 2.0 * PI) / 3.0) - a3);

    // StackVector begin() and end() are pointers, so std::sort works directly.
    std::sort(res.begin(), res.end());
    return res;
  } else {
    double inner = std::abs(R) + std::sqrt(R2 - Q3);
    double A = (R < 0 ? 1.0 : -1.0) * std::cbrt(inner);
    double B = std::abs(A) >= epsilon ? Q / A : 0.0;
    res.push_back(A + B - a3);
    return res;
  }
}

Roots3 GeometryEpsilon::solveCubic(double a, double b, double c,
                                   double d) const {
  if (std::abs(a) < epsilon) {
    Roots3 res;
    if (std::abs(b) < epsilon) {
      if (std::abs(c) < epsilon) {
        if (std::abs(d) < epsilon)
          res.push_back(0.0);
        return res;
      }
      res.push_back(-d / c);
      return res;
    }
    double b2 = 2.0 * b;
    double D = c * c - 4.0 * b * d;
    if (std::abs(D) < epsilon) {
      res.push_back(-c / b2);
      return res;
    } else if (D > 0) {
      D = std::sqrt(D);
      res.push_back((-c + D) / b2);
      res.push_back((-c - D) / b2);
      std::sort(res.begin(), res.end());
      return res;
    }
    return res;
  }
  return solveCubicNormalized(b / a, c / a, d / a);
}

// -----------------------------------------------------------------
// Low-level geometry implementation for segments.
// -----------------------------------------------------------------

StackVector<double, 8> SegmentCurve::boundingTValues() const {
  SegmentTValuesBuilder result(geo);
  auto bounds = [this, &result](double x0, double x1, double x2, double x3) {
    double a = 3.0 * x3 - 9.0 * x2 + 9.0 * x1 - 3.0 * x0;
    double b = 6.0 * x0 - 12.0 * x1 + 6.0 * x2;
    double c = 3.0 * x1 - 3.0 * x0;
    if (geo->snap0(a) == 0.0) {
      if (geo->snap0(b) != 0.0)
        result.add(-c / b);
    } else {
      double disc = b * b - 4.0 * a * c;
      if (disc >= 0.0) {
        double sq = std::sqrt(disc);
        result.add((-b + sq) / (2.0 * a));
        result.add((-b - sq) / (2.0 * a));
      }
    }
  };

  bounds(p0.x, p1.x, p2.x, p3.x);
  bounds(p0.y, p1.y, p2.y, p3.y);
  return result.list();
}

StackVector<double, 8> SegmentCurve::inflectionTValues() const {
  SegmentTValuesBuilder result(geo);
  result.addArray(boundingTValues());

  double p10x = 3.0 * (p1.x - p0.x);
  double p10y = 3.0 * (p1.y - p0.y);
  double p21x = 6.0 * (p2.x - p1.x);
  double p21y = 6.0 * (p2.y - p1.y);
  double p32x = 3.0 * (p3.x - p2.x);
  double p32y = 3.0 * (p3.y - p2.y);
  double p210x = 6.0 * (p2.x - 2.0 * p1.x + p0.x);
  double p210y = 6.0 * (p2.y - 2.0 * p1.y + p0.y);
  double p321x = 6.0 * (p3.x - 2.0 * p2.x + p1.x);
  double p321y = 6.0 * (p3.y - 2.0 * p2.y + p1.y);

  double qx = p10x - p21x + p32x;
  double qy = p10y - p21y + p32y;
  double rx = p21x - 2.0 * p10x;
  double ry = p21y - 2.0 * p10y;
  double sx = p10x;
  double sy = p10y;
  double ux = p321x - p210x;
  double uy = p321y - p210y;
  double vx = p210x;
  double vy = p210y;

  double A = qx * uy - qy * ux;
  double B = qx * vy + rx * uy - qy * vx - ry * ux;
  double C = rx * vy + sx * uy - ry * vx - sy * ux;
  double D = sx * vy - sy * vx;

  for (double s : geo->solveCubic(A, B, C, D)) {
    result.add(s);
  }
  return result.list();
}

BoundingBox SegmentCurve::boundingBox() const {
  Vec2 minVal = {std::min(p0.x, p3.x), std::min(p0.y, p3.y)};
  Vec2 maxVal = {std::max(p0.x, p3.x), std::max(p0.y, p3.y)};
  for (double t : boundingTValues()) {
    Vec2 p = point(t);
    minVal.x = std::min(minVal.x, p.x);
    minVal.y = std::min(minVal.y, p.y);
    maxVal.x = std::max(maxVal.x, p.x);
    maxVal.y = std::max(maxVal.y, p.y);
  }
  return {minVal, maxVal};
}

std::optional<double> SegmentCurve::mapXtoT(double x, bool force) const {
  if (geo->snap0(p0.x - x) == 0.0)
    return 0.0;
  if (geo->snap0(p3.x - x) == 0.0)
    return 1.0;

  double px0 = p0.x - x;
  double px1 = p1.x - x;
  double px2 = p2.x - x;
  double px3 = p3.x - x;

  std::vector<double> R = {px3 - 3.0 * px2 + 3.0 * px1 - px0,
                           3.0 * px2 - 6.0 * px1 + 3.0 * px0,
                           3.0 * px1 - 3.0 * px0, px0};

  for (double t : geo->solveCubic(R[0], R[1], R[2], R[3])) {
    double ts = geo->snap01(t);
    if (ts >= 0.0 && ts <= 1.0)
      return t;
  }

  if (force || (x >= std::min(p0.x, p3.x) && x <= std::max(p0.x, p3.x))) {
    for (int attempt = 0; attempt < 4; attempt++) {
      int ii = -1;
      for (int i = 0; i < 4; i++) {
        if (R[i] != 0.0 && (ii < 0 || std::abs(R[i]) < std::abs(R[ii]))) {
          ii = i;
        }
      }
      if (ii < 0)
        return 0.0;
      R[ii] = 0.0;

      for (double t : geo->solveCubic(R[0], R[1], R[2], R[3])) {
        double ts = geo->snap01(t);
        if (ts >= 0.0 && ts <= 1.0)
          return t;
      }
    }
  }
  return std::nullopt;
}

std::optional<double> SegmentCurve::mapXtoY(double x, bool force) const {
  auto t = mapXtoT(x, force);
  if (!t.has_value())
    return std::nullopt;
  return point(t.value()).y;
}

bool SegmentCurve::pointOn(const Vec2 &p) const {
  if (geo->isEqualVec2(p0, p) || geo->isEqualVec2(p3, p))
    return true;
  auto y = mapXtoY(p.x);
  if (!y.has_value())
    return false;
  return geo->snap0(y.value() - p.y) == 0.0;
}

std::optional<SegmentLine> SegmentCurve::toLine() const {
  if ((geo->snap0(p0.x - p1.x) == 0.0 && geo->snap0(p0.x - p2.x) == 0.0 &&
       geo->snap0(p0.x - p3.x) == 0.0) ||
      (geo->snap0(p0.y - p1.y) == 0.0 && geo->snap0(p0.y - p2.y) == 0.0 &&
       geo->snap0(p0.y - p3.y) == 0.0)) {
    return SegmentLine(p0, p3, geo);
  }
  return std::nullopt;
}

double projectPointOntoSegmentLine(const Vec2 &p, const SegmentLine &seg) {
  double dx = seg.p1.x - seg.p0.x;
  double dy = seg.p1.y - seg.p0.y;
  double px = p.x - seg.p0.x;
  double py = p.y - seg.p0.y;
  double dist = dx * dx + dy * dy;
  double dot = px * dx + py * dy;
  return dot / dist;
}

IntersectionResult segmentLineIntersectSegmentLine(const SegmentLine &segA,
                                                   const SegmentLine &segB,
                                                   bool allowOutOfRange) {
  const Geometry *geo = segA.geo;
  Vec2 a0 = segA.p0;
  Vec2 a1 = segA.p1;
  Vec2 b0 = segB.p0;
  Vec2 b1 = segB.p1;

  double adx = a1.x - a0.x;
  double ady = a1.y - a0.y;
  double bdx = b1.x - b0.x;
  double bdy = b1.y - b0.y;

  double axb = adx * bdy - ady * bdx;
  if (geo->snap0(axb) == 0.0) {
    if (!geo->isCollinear(a0, a1, b0)) {
      return std::monostate{};
    }

    double tB0onA = projectPointOntoSegmentLine(segB.p0, segA);
    double tB1onA = projectPointOntoSegmentLine(segB.p1, segA);
    double tAMin = geo->snap01(std::min(tB0onA, tB1onA));
    double tAMax = geo->snap01(std::max(tB0onA, tB1onA));
    if (tAMax < 0.0 || tAMin > 1.0)
      return std::monostate{};

    double tA0onB = projectPointOntoSegmentLine(segA.p0, segB);
    double tA1onB = projectPointOntoSegmentLine(segA.p1, segB);
    double tBMin = geo->snap01(std::min(tA0onB, tA1onB));
    double tBMax = geo->snap01(std::max(tA0onB, tA1onB));
    if (tBMax < 0.0 || tBMin > 1.0)
      return std::monostate{};

    return SegmentTRangePairs{{std::max(0.0, tAMin), std::max(0.0, tBMin)},
                              {std::min(1.0, tAMax), std::min(1.0, tBMax)}};
  }

  double dx = a0.x - b0.x;
  double dy = a0.y - b0.y;
  return SegmentTValuePairsBuilder(allowOutOfRange, geo)
      .add((bdx * dy - bdy * dx) / axb, (adx * dy - ady * dx) / axb)
      .done();
}

IntersectionResult segmentLineIntersectSegmentCurve(const SegmentLine &segA,
                                                    const SegmentCurve &segB,
                                                    bool allowOutOfRange,
                                                    bool invert) {
  const Geometry *geo = segA.geo;
  Vec2 a0 = segA.p0;
  Vec2 a1 = segA.p1;

  double A = a1.y - a0.y;
  double B = a0.x - a1.x;

  if (geo->snap0(B) == 0.0) {
    auto t_opt = segB.mapXtoT(a0.x, false);
    if (!t_opt.has_value())
      return std::monostate{};
    double t = t_opt.value();
    double y = segB.point(t).y;
    double s = (y - a0.y) / A;
    SegmentTValuePairsBuilder result(allowOutOfRange, geo);
    if (invert)
      result.add(t, s);
    else
      result.add(s, t);
    return result.done();
  }

  double C = A * a0.x + B * a0.y;
  std::array<double, 4> bx = segB.getCubicCoefficients(0);
  std::array<double, 4> by = segB.getCubicCoefficients(1);

  double rA = A * bx[0] + B * by[0];
  double rB = A * bx[1] + B * by[1];
  double rC = A * bx[2] + B * by[2];
  double rD = A * bx[3] + B * by[3] - C;

  Roots3 roots = geo->solveCubic(rA, rB, rC, rD);
  SegmentTValuePairsBuilder result(allowOutOfRange, geo);

  if (geo->snap0(A) == 0.0) {
    for (double t : roots) {
      double X = bx[0] * t * t * t + bx[1] * t * t + bx[2] * t + bx[3];
      double s = (a0.x - X) / B;
      if (invert)
        result.add(t, s);
      else
        result.add(s, t);
    }
  } else {
    for (double t : roots) {
      double Y = by[0] * t * t * t + by[1] * t * t + by[2] * t + by[3];
      double s = (Y - a0.y) / A;
      if (invert)
        result.add(t, s);
      else
        result.add(s, t);
    }
  }

  return result.done();
}

IntersectionResult segmentCurveIntersectSegmentCurve(const SegmentCurve &segA,
                                                     const SegmentCurve &segB,
                                                     bool allowOutOfRange) {
  const Geometry *geo = segA.geo;

  if (geo->isEqualVec2(segA.p0, segB.p0)) {
    if (geo->isEqualVec2(segA.p3, segB.p3)) {
      if (geo->isEqualVec2(segA.p1, segB.p1) &&
          geo->isEqualVec2(segA.p2, segB.p2)) {
        return SegmentTRangePairs{{0.0, 0.0}, {1.0, 1.0}};
      } else {
        SegmentTValuePairs res;
        res.tValuePairs.push_back({0.0, 0.0});
        res.tValuePairs.push_back({1.0, 1.0});
        return res;
      }
    } else {
      SegmentTValuePairs res;
      res.tValuePairs.push_back({0.0, 0.0});
      return res;
    }
  } else if (geo->isEqualVec2(segA.p0, segB.p3)) {
    SegmentTValuePairs res;
    res.tValuePairs.push_back({0.0, 1.0});
    return res;
  } else if (geo->isEqualVec2(segA.p3, segB.p0)) {
    SegmentTValuePairs res;
    res.tValuePairs.push_back({1.0, 0.0});
    return res;
  } else if (geo->isEqualVec2(segA.p3, segB.p3)) {
    SegmentTValuePairs res;
    res.tValuePairs.push_back({1.0, 1.0});
    return res;
  }

  SegmentTValuePairsBuilder result(allowOutOfRange, geo);

  // Y-combinator style recursive lambda for bounding box checks
  auto checkCurves = [&](auto &self, const SegmentCurve &c1, double t1L,
                         double t1R, const SegmentCurve &c2, double t2L,
                         double t2R) -> void {
    BoundingBox bbox1 = c1.boundingBox();
    BoundingBox bbox2 = c2.boundingBox();

    if (!boundingBoxesIntersect(bbox1, bbox2))
      return;

    double t1M = (t1L + t1R) / 2.0;
    double t2M = (t2L + t2R) / 2.0;

    if (geo->snap0(t1R - t1L) == 0.0 && geo->snap0(t2R - t2L) == 0.0) {
      result.add(t1M, t2M);
      return;
    }
    StackVector<SegmentCurve, 2> splitC1;
    StackVector<SegmentCurve, 2> splitC2;
    c1.splitSingle(0.5, splitC1);
    c2.splitSingle(0.5, splitC2);

    // Guard against infinite loops if split falls back to the original segment.
    assert(splitC1.size() == 2 && splitC2.size() == 2 &&
           "Split must return exactly 2 curves");

    self(self, splitC1[0], t1L, t1M, splitC2[0], t2L, t2M);
    self(self, splitC1[1], t1M, t1R, splitC2[0], t2L, t2M);
    self(self, splitC1[0], t1L, t1M, splitC2[1], t2M, t2R);
    self(self, splitC1[1], t1M, t1R, splitC2[1], t2M, t2R);
  };

  checkCurves(checkCurves, segA, 0.0, 1.0, segB, 0.0, 1.0);
  return result.done();
}

IntersectionResult segmentsIntersect(const Segment &segA, const Segment &segB,
                                     bool allowOutOfRange) {
  return std::visit(
      [&](auto &&argA, auto &&argB) -> IntersectionResult {
        using T = std::decay_t<decltype(argA)>;
        using U = std::decay_t<decltype(argB)>;

        if constexpr (std::is_same_v<T, SegmentLine> &&
                      std::is_same_v<U, SegmentLine>) {
          return segmentLineIntersectSegmentLine(argA, argB, allowOutOfRange);
        } else if constexpr (std::is_same_v<T, SegmentLine> &&
                             std::is_same_v<U, SegmentCurve>) {
          return segmentLineIntersectSegmentCurve(argA, argB, allowOutOfRange,
                                                  false);
        } else if constexpr (std::is_same_v<T, SegmentCurve> &&
                             std::is_same_v<U, SegmentLine>) {
          return segmentLineIntersectSegmentCurve(argB, argA, allowOutOfRange,
                                                  true);
        } else if constexpr (std::is_same_v<T, SegmentCurve> &&
                             std::is_same_v<U, SegmentCurve>) {
          return segmentCurveIntersectSegmentCurve(argA, argB, allowOutOfRange);
        } else {
          throw std::runtime_error(
              "PolyBool: Unknown segment instance in segmentsIntersect");
        }
      },
      segA, segB);
}

// --- Arena Factory Method Implementations ---
SegmentBool *Intersecter::makeSegmentBool(const Segment &data,
                                          const SegmentBoolFill &fill,
                                          bool closed) {
  segmentArena.emplace_back(data, fill, closed);
  return &segmentArena.back();
}

EventBool *Intersecter::makeEventBool(bool isStart, const Vec2 &p,
                                      SegmentBool *seg, bool primary) {
  eventArena.emplace_back(isStart, p, seg, primary);
  return &eventArena.back();
}

// ============================================================================
// Anonymous namespace: helpers for accessing std::variant<SegmentLine, SegmentCurve>.
// ============================================================================
namespace {
auto getStart = [](const Segment &s) -> Vec2 {
  return std::visit([](auto &&arg) { return arg.start(); }, s);
};
auto getStart2 = [](const Segment &s) -> Vec2 {
  return std::visit([](auto &&arg) { return arg.start2(); }, s);
};
auto getEnd = [](const Segment &s) -> Vec2 {
  return std::visit([](auto &&arg) { return arg.end(); }, s);
};
auto setStart = [](Segment &s, const Vec2 &p) {
  std::visit([&p](auto &&arg) { arg.setStart(p); }, s);
};
auto setEnd = [](Segment &s, const Vec2 &p) {
  std::visit([&p](auto &&arg) { arg.setEnd(p); }, s);
};
auto pointOn = [](const Segment &s, const Vec2 &p) -> bool {
  return std::visit([&p](auto &&arg) { return arg.pointOn(p); }, s);
};
auto getPoint = [](const Segment &s, double t) -> Vec2 {
  return std::visit([t](auto &&arg) { return arg.point(t); }, s);
};

// Simple sign function equivalent to JavaScript Math.sign.
auto math_sign = [](double val) -> int { return (val > 0.0) - (val < 0.0); };
} // namespace

void Intersecter::beginPath() { currentPath.clear(); }

void Intersecter::closePath() {
  for (SegmentBool *seg : currentPath) {
    seg->closed = true;
  }
}

EventBool *Intersecter::addSegment(SegmentBool *seg, bool primary) {
  EventBool *evStart = makeEventBool(true, getStart(seg->data), seg, primary);
  EventBool *evEnd = makeEventBool(false, getEnd(seg->data), seg, primary);

  evStart->other = evEnd;
  evEnd->other = evStart;

  addEvent(evStart);
  addEvent(evEnd);
  return evStart;
}

void Intersecter::addLine(const Vec2 &from, const Vec2 &to, bool primary) {
  int f = geo->compareVec2(from, to);
  if (f == 0)
    return; // Skip zero-length line segments.

  Segment segData = SegmentLine(f < 0 ? from : to, f < 0 ? to : from, geo);
  SegmentBool *seg = makeSegmentBool(segData, {}, false);
  currentPath.push_back(seg);
  addSegment(seg, primary);
}

void Intersecter::addCurve(const Vec2 &from, const Vec2 &c1, const Vec2 &c2,
                           const Vec2 &to, bool primary) {
  SegmentCurve original(from, c1, c2, to, geo);
  std::vector<SegmentCurve> curves;
  auto infs = original.inflectionTValues();
  original.split(infs.data, infs.count, curves);

  for (const auto &curve : curves) {
    int f = geo->compareVec2(curve.start(), curve.end());
    if (f == 0)
      continue;

    auto lineOpt = curve.toLine();
    if (lineOpt.has_value()) {
      addLine(lineOpt->p0, lineOpt->p1, primary);
    } else {
      Segment segData = f < 0 ? curve : curve.reverse();
      SegmentBool *seg = makeSegmentBool(segData, {}, false);
      currentPath.push_back(seg);
      addSegment(seg, primary);
    }
  }
}

void Intersecter::addEvent(EventBool *ev) {
  events.insertBefore(ev, [this, ev](EventBool *here) -> int {
    if (here == ev)
      return 0;
    return this->compareEvents(ev->isStart, ev->p, ev->other->p, ev->seg->data,
                               here->isStart, here->p, here->other->p,
                               here->seg->data);
  });
}

int Intersecter::compareEvents(bool aStart, const Vec2 &a1, const Vec2 &a2,
                               const Segment &aSeg, bool bStart, const Vec2 &b1,
                               const Vec2 &b2, const Segment &bSeg) {
  int comp = geo->compareVec2(a1, b1);
  if (comp != 0)
    return comp;

  bool aIsLine = std::holds_alternative<SegmentLine>(aSeg);
  bool bIsLine = std::holds_alternative<SegmentLine>(bSeg);

  if (aIsLine && bIsLine && geo->isEqualVec2(a2, b2)) {
    return 0;
  }

  if (aStart != bStart) {
    return aStart ? 1 : -1;
  }

  return compareSegments(bSeg, aSeg);
}

EventBool *Intersecter::divideEvent(EventBool *ev, double t, const Vec2 &p) {
  splitBuffer.clear();
  // 1. Split the segment.
  std::visit([&](auto &&arg) { arg.splitSingle(t, splitBuffer); },
             ev->seg->data);
  // Safety checks.
  if (splitBuffer.size() < 2) {
    throw std::runtime_error(
        "PolyBool: Split failed to produce at least 2 segments.");
  }
  Segment left = splitBuffer[0];
  Segment right = splitBuffer[1];

  // 2. Force the exact intersection point.
  setEnd(left, p);
  setStart(right, p);

  // 3. Create and register the new object.
  SegmentBool *ns = makeSegmentBool(right, ev->seg->myFill, ev->seg->closed);

  // 4. Update the original object in the arena.
  events.remove(ev->other);
  ev->seg->data = left;
  ev->other->p = p;
  addEvent(ev->other);

  return addSegment(ns, ev->primary);
}

int Intersecter::compareSegments(const Segment &seg1, const Segment &seg2) {
  Vec2 A = getStart(seg1);
  Vec2 B = getStart2(seg2);
  Vec2 C = getStart(seg2);

  if (pointOn(seg2, A)) {
    A = getStart2(seg1);
    if (pointOn(seg2, A)) {
      bool s1Line = std::holds_alternative<SegmentLine>(seg1);
      bool s2Line = std::holds_alternative<SegmentLine>(seg2);
      bool s1Curve = std::holds_alternative<SegmentCurve>(seg1);
      bool s2Curve = std::holds_alternative<SegmentCurve>(seg2);

      if (s1Line) {
        if (s2Line)
          return 0;
        if (s2Curve)
          A = std::get<SegmentLine>(seg1).point(0.5);
      }
      if (s1Curve)
        A = getEnd(seg1);
    }
    if (std::holds_alternative<SegmentCurve>(seg2)) {
      if (geo->snap0(A.x - C.x) == 0.0 && geo->snap0(B.x - C.x) == 0.0) {
        return math_sign(C.y - A.y);
      }
    }
  } else {
    if (std::holds_alternative<SegmentCurve>(seg2)) {
      auto optY = std::get<SegmentCurve>(seg2).mapXtoY(A.x, true);
      if (optY.has_value())
        return math_sign(optY.value() - A.y);
    }
    if (std::holds_alternative<SegmentCurve>(seg1)) {
      auto i = segmentsIntersect(seg1, seg2, true);
      if (std::holds_alternative<SegmentTValuePairs>(i)) {
        for (const auto &pair : std::get<SegmentTValuePairs>(i).tValuePairs) {
          double t = geo->snap01(pair.tA);
          if (t > 0.0 && t < 1.0) {
            B = std::get<SegmentCurve>(seg1).point(t);
            break;
          }
        }
      }
    }
  }

  double cross = (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);
  return math_sign(cross);
}

ListBoolTransition<EventBool *>
Intersecter::statusFindSurrounding(EventBool *ev) {
  return status.findTransition(ev, [this, ev](EventBool *here) -> int {
    if (ev == here)
      return 0;
    int c = this->compareSegments(ev->seg->data, here->seg->data);
    return c == 0 ? -1 : c;
  });
}

EventBool *Intersecter::checkIntersection(EventBool *ev1, EventBool *ev2) {
  SegmentBool *seg1 = ev1->seg;
  SegmentBool *seg2 = ev2->seg;

  IntersectionResult i = segmentsIntersect(seg1->data, seg2->data, false);

  if (std::holds_alternative<std::monostate>(i)) {
    return nullptr;
  } else if (std::holds_alternative<SegmentTRangePairs>(i)) {
    auto range = std::get<SegmentTRangePairs>(i);
    double tA1 = range.tStart.tA, tB1 = range.tStart.tB;
    double tA2 = range.tEnd.tA, tB2 = range.tEnd.tB;

    if ((tA1 == 1.0 && tA2 == 1.0 && tB1 == 0.0 && tB2 == 0.0) ||
        (tA1 == 0.0 && tA2 == 0.0 && tB1 == 1.0 && tB2 == 1.0)) {
      return nullptr;
    }

    if (tA1 == 0.0 && tA2 == 1.0 && tB1 == 0.0 && tB2 == 1.0) {
      return ev2;
    }

    Vec2 a1 = getStart(seg1->data);
    Vec2 a2 = getEnd(seg1->data);
    Vec2 b2 = getEnd(seg2->data);

    if (tA1 == 0.0 && tB1 == 0.0) {
      if (tA2 == 1.0)
        divideEvent(ev2, tB2, a2);
      else
        divideEvent(ev1, tA2, b2);
      return ev2;
    } else if (tB1 > 0.0 && tB1 < 1.0) {
      if (tA2 == 1.0 && tB2 == 1.0) {
        divideEvent(ev2, tB1, a1);
      } else {
        if (tA2 == 1.0)
          divideEvent(ev2, tB2, a2);
        else
          divideEvent(ev1, tA2, b2);
        divideEvent(ev2, tB1, a1);
      }
    }
    return nullptr;

  } else if (std::holds_alternative<SegmentTValuePairs>(i)) {
    auto pairs = std::get<SegmentTValuePairs>(i).tValuePairs;
    if (pairs.empty())
      return nullptr;

    TValuePair minPair = pairs[0];
    for (size_t j = 1; j < pairs.size(); ++j) {
      bool skip = (minPair.tA == 0.0 && minPair.tB == 0.0) ||
                  (minPair.tA == 0.0 && minPair.tB == 1.0) ||
                  (minPair.tA == 1.0 && minPair.tB == 0.0) ||
                  (minPair.tA == 1.0 && minPair.tB == 1.0);
      if (!skip)
        break;
      minPair = pairs[j];
    }

    double tA = minPair.tA;
    double tB = minPair.tB;

    Vec2 p;
    if (tB == 0.0)
      p = getStart(seg2->data);
    else if (tB == 1.0)
      p = getEnd(seg2->data);
    else if (tA == 0.0)
      p = getStart(seg1->data);
    else if (tA == 1.0)
      p = getEnd(seg1->data);
    else
      p = getPoint(seg1->data, tA);

    if (tA > 0.0 && tA < 1.0)
      divideEvent(ev1, tA, p);
    if (tB > 0.0 && tB < 1.0)
      divideEvent(ev2, tB, p);

    return nullptr;
  }

  throw std::runtime_error("PolyBool: Unknown intersection type");
}

// --- Core calculate() Method ---
std::vector<SegmentBool> Intersecter::calculate() {
  std::vector<SegmentBool> resultSegments;

  while (!events.isEmpty()) {
    EventBool *ev = events.getHead();

    if (ev->isStart) {
      auto surrounding = statusFindSurrounding(ev);
      EventBool *above = surrounding.before;
      EventBool *below = surrounding.after;

      auto checkBothIntersections = [&]() -> EventBool * {
        if (above) {
          EventBool *eve = checkIntersection(ev, above);
          if (eve)
            return eve;
        }
        if (below)
          return checkIntersection(ev, below);
        return nullptr;
      };

      EventBool *eve = checkBothIntersections();
      if (eve) {
        // Overlap merge logic.
        if (selfIntersection) {
          bool toggle = (ev->seg->myFill.below == std::nullopt)
                            ? ev->seg->closed
                            : (ev->seg->myFill.above != ev->seg->myFill.below);
          if (toggle)
            eve->seg->myFill.above = !eve->seg->myFill.above.value_or(false);
        } else {
          eve->seg->otherFill = ev->seg->myFill;
        }
        events.remove(ev->other);
        events.remove(ev);
      }

      if (events.getHead() != ev)
        continue;

      // Calculate fill flags.
      if (selfIntersection) {
        bool toggle;
        if (!ev->seg->myFill.below.has_value()) {
          toggle = ev->seg->closed;
        } else {
          toggle = ev->seg->myFill.above.value_or(false) !=
                   ev->seg->myFill.below.value_or(false);
        }

        if (!below) {
          ev->seg->myFill.below = false;
        } else {
          ev->seg->myFill.below = below->seg->myFill.above;
        }

        ev->seg->myFill.above = toggle ? !ev->seg->myFill.below.value()
                                       : ev->seg->myFill.below.value();

      } else {

        if (!ev->seg->otherFill.has_value()) {
          bool inside = false;
          if (!below) {
            inside = false;
          } else {
            if (ev->primary == below->primary) {
              if (!below->seg->otherFill.has_value()) {
                throw std::runtime_error(
                    "PolyBool: Unexpected state of otherFill (null)");
              }
              inside = below->seg->otherFill->above.value_or(false);
            } else {
              inside = below->seg->myFill.above.value_or(false);
            }
          }
          ev->seg->otherFill = SegmentBoolFill{inside, inside};
        }
      }

      ev->other->status = surrounding.insert(ev);

    } else {
      // End-event handling.
      EventBool *st = ev->status;
      if (!st)
        throw std::runtime_error("PolyBool: Zero-length segment detected; your "
                                 "epsilon is probably too small or too large");

      int i = status.getIndex(st);
      if (i > 0 && i < (int)status.nodes.size() - 1) {
        EventBool *before = status.nodes[i - 1];
        EventBool *after = status.nodes[i + 1];
        checkIntersection(before, after);
      }

      status.remove(st);

      if (!ev->primary) {
        auto s = ev->seg->myFill;
        ev->seg->myFill = ev->seg->otherFill.value_or(SegmentBoolFill{});
        ev->seg->otherFill = s;
      }

      resultSegments.push_back(*(ev->seg));
    }

    events.removeHead();
  }

  return resultSegments;
}

// --- Merge Helper Functions ---

std::optional<SegmentLine> joinLines(const SegmentLine &seg1,
                                     const SegmentLine &seg2,
                                     const Geometry *geo) {
  if (geo->isCollinear(seg1.p0, seg1.p1, seg2.p1)) {
    return SegmentLine(seg1.p0, seg2.p1, geo);
  }
  return std::nullopt;
}

std::optional<SegmentCurve> joinCurves(const SegmentCurve &seg1,
                                       const SegmentCurve &seg2,
                                       const Geometry *geo) {
  if (geo->isCollinear(seg1.p2, seg1.p3, seg2.p1)) {
    double dx = seg2.p1.x - seg1.p2.x;
    double dy = seg2.p1.y - seg1.p2.y;
    double t = std::abs(dx) > std::abs(dy) ? (seg1.p3.x - seg1.p2.x) / dx
                                           : (seg1.p3.y - seg1.p2.y) / dy;

    double ts = geo->snap01(t);
    if (ts != 0.0 && ts != 1.0) {
      SegmentCurve ns(seg1.p0,
                      {seg1.p0.x + (seg1.p1.x - seg1.p0.x) / t,
                       seg1.p0.y + (seg1.p1.y - seg1.p0.y) / t},
                      {seg2.p2.x - (t * (seg2.p3.x - seg2.p2.x)) / (1.0 - t),
                       seg2.p2.y - (t * (seg2.p3.y - seg2.p2.y)) / (1.0 - t)},
                      seg2.p3, geo);

      std::vector<SegmentCurve> split_result;
      ns.splitSingle(t, split_result);
      if (split_result.size() == 2 && split_result[0].isEqual(seg1) &&
          split_result[1].isEqual(seg2)) {
        return ns;
      }
    }
  }
  return std::nullopt;
}

std::optional<Segment> joinSegments(const Segment &seg1, const Segment &seg2,
                                    const Geometry *geo) {
  if (&seg1 == &seg2)
    return std::nullopt;
  // std::visit provides pattern matching over variants.
  return std::visit(
      [&](auto &&arg1, auto &&arg2) -> std::optional<Segment> {
        using T1 = std::decay_t<decltype(arg1)>;
        using T2 = std::decay_t<decltype(arg2)>;

        if constexpr (std::is_same_v<T1, SegmentLine> &&
                      std::is_same_v<T2, SegmentLine>) {
          auto res = joinLines(arg1, arg2, geo);
          if (res)
            return *res;
        } else if constexpr (std::is_same_v<T1, SegmentCurve> &&
                             std::is_same_v<T2, SegmentCurve>) {
          auto res = joinCurves(arg1, arg2, geo);
          if (res)
            return *res;
        }
        return std::nullopt;
      },
      seg1, seg2);
}

// --- Internal Chain Storage ---
struct SegsFill {
  std::deque<Segment> segs; // deque supports efficient front/back insertion.
  bool fill;
};

// --- Main Chainer Function ---

std::vector<std::vector<Segment>>
SegmentChainer(const std::vector<SegmentBool> &segments, const Geometry *geo) {
  std::vector<SegsFill> closedChains;
  std::vector<SegsFill> openChains;
  std::vector<std::vector<Segment>> regions;

  auto reverseSeg = [](const Segment &s) -> Segment {
    return std::visit([](auto &&arg) -> Segment { return arg.reverse(); }, s);
  };

  auto getStart = [](const Segment &s) -> Vec2 {
    return std::visit([](auto &&arg) { return arg.start(); }, s);
  };

  auto getEnd = [](const Segment &s) -> Vec2 {
    return std::visit([](auto &&arg) { return arg.end(); }, s);
  };

  for (const auto &segb : segments) {
    Segment seg = segb.data;
    bool is_closed = segb.closed;
    auto &chains = is_closed ? closedChains : openChains;
    Vec2 pt1 = getStart(seg);
    Vec2 pt2 = getEnd(seg);

    auto reverseChain = [&](size_t index) -> std::deque<Segment> & {
      std::deque<Segment> newChain;
      for (const auto &s : chains[index].segs) {
        newChain.push_front(reverseSeg(s));
      }
      chains[index].segs = std::move(newChain);
      chains[index].fill = !chains[index].fill;
      return chains[index].segs;
    };

    // Skip zero-length line segments.
    if (std::holds_alternative<SegmentLine>(seg)) {
      if (geo->isEqualVec2(pt1, pt2)) {
        std::cerr << "PolyBool: Warning: Zero-length segment detected\n";
        continue;
      }
    }
    struct Match {
      size_t index = 0;
      bool matchesHead = false;
      bool matchesPt1 = false;
    };

    Match firstMatch, secondMatch;
    Match *nextMatch = &firstMatch;

    auto setMatch = [&](size_t index, bool matchesHead, bool matchesPt1) {
      if (nextMatch) {
        nextMatch->index = index;
        nextMatch->matchesHead = matchesHead;
        nextMatch->matchesPt1 = matchesPt1;
      }
      if (nextMatch == &firstMatch) {
        nextMatch = &secondMatch;
        return false;
      }
      nextMatch = nullptr;
      return true;
    };

    // Find connected chains.
    for (size_t i = 0; i < chains.size(); ++i) {
      const auto &chain = chains[i].segs;
      Vec2 head = getStart(chain.front());
      Vec2 tail = getEnd(chain.back());

      if (geo->isEqualVec2(head, pt1)) {
        if (setMatch(i, true, true))
          break;
      } else if (geo->isEqualVec2(head, pt2)) {
        if (setMatch(i, true, false))
          break;
      } else if (geo->isEqualVec2(tail, pt1)) {
        if (setMatch(i, false, true))
          break;
      } else if (geo->isEqualVec2(tail, pt2)) {
        if (setMatch(i, false, false))
          break;
      }
    }

    if (nextMatch == &firstMatch) {
      // No match; create a new chain.
      chains.push_back(
          {std::deque<Segment>{seg}, segb.myFill.above.value_or(false)});
    } else if (nextMatch == &secondMatch) {
      // Matched one end of a chain.
      size_t index = firstMatch.index;
      auto &chain_data = chains[index];
      auto &chain = chain_data.segs;

      if (firstMatch.matchesHead) {
        if (firstMatch.matchesPt1)
          seg = reverseSeg(seg);
        chain.push_front(seg);

        // Simplify the head.
        if (chain.size() >= 2) {
          auto newSeg = joinSegments(chain[0], chain[1], geo);
          if (newSeg) {
            chain.pop_front();
            chain[0] = *newSeg;
          }
        }
      } else {
        if (!firstMatch.matchesPt1)
          seg = reverseSeg(seg);
        chain.push_back(seg);

        // Simplify the tail.
        if (chain.size() >= 2) {
          auto newSeg =
              joinSegments(chain[chain.size() - 2], chain.back(), geo);
          if (newSeg) {
            chain.pop_back();
            chain.back() = *newSeg;
          }
        }
      }

      // Check for closure.
      if (is_closed && chain.size() > 0 &&
          geo->isEqualVec2(getStart(chain.front()), getEnd(chain.back()))) {
        double winding = 0.0;
        Vec2 last = getStart(chain.front());
        for (const auto &s : chain) {
          Vec2 here = getEnd(s);
          winding += here.y * last.x - here.x * last.y;
          last = here;
        }
        bool isClockwise = winding < 0.0;
        if (isClockwise == chain_data.fill) {
          reverseChain(index);
        }

        auto newStart = joinSegments(chain.back(), chain.front(), geo);
        if (newStart) {
          chain.pop_back();
          chain.front() = *newStart;
        }

        regions.push_back(std::vector<Segment>(chain.begin(), chain.end()));
        chains.erase(chains.begin() + index);
      }

    } else {
      // Matched two chains; bridge them.
      auto appendChain = [&](size_t idx1, size_t idx2, Segment s) {
        auto &chain1 = chains[idx1].segs;
        auto &chain2 = chains[idx2].segs;

        chain1.push_back(s);
        auto newEnd =
            joinSegments(chain1[chain1.size() - 2], chain1.back(), geo);
        if (newEnd) {
          chain1.pop_back();
          chain1.back() = *newEnd;
        }

        auto newJoin = joinSegments(chain1.back(), chain2.front(), geo);
        if (newJoin) {
          chain2.pop_front();
          chain1.back() = *newJoin;
        }

        chain1.insert(chain1.end(), std::make_move_iterator(chain2.begin()),
                      std::make_move_iterator(chain2.end()));
        // Mark for deletion and process later.
        chains[idx2].segs.clear();
      };

      size_t F = firstMatch.index;
      size_t S = secondMatch.index;
      bool reverseF = chains[F].segs.size() < chains[S].segs.size();

      if (firstMatch.matchesHead) {
        if (secondMatch.matchesHead) {
          if (reverseF) {
            if (!firstMatch.matchesPt1)
              seg = reverseSeg(seg);
            reverseChain(F);
            appendChain(F, S, seg);
          } else {
            if (firstMatch.matchesPt1)
              seg = reverseSeg(seg);
            reverseChain(S);
            appendChain(S, F, seg);
          }
        } else {
          if (firstMatch.matchesPt1)
            seg = reverseSeg(seg);
          appendChain(S, F, seg);
        }
      } else {
        if (secondMatch.matchesHead) {
          if (!firstMatch.matchesPt1)
            seg = reverseSeg(seg);
          appendChain(F, S, seg);
        } else {
          if (reverseF) {
            if (firstMatch.matchesPt1)
              seg = reverseSeg(seg);
            reverseChain(F);
            appendChain(S, F, seg);
          } else {
            if (!firstMatch.matchesPt1)
              seg = reverseSeg(seg);
            reverseChain(S);
            appendChain(F, S, seg);
          }
        }
      }

      // Remove empty chains that were merged away.
      auto it =
          std::remove_if(chains.begin(), chains.end(),
                         [](const SegsFill &sf) { return sf.segs.empty(); });
      chains.erase(it, chains.end());
    }
  }

  // Add remaining open chains to the result.
  for (const auto &chain : openChains) {
    regions.push_back(
        std::vector<Segment>(chain.segs.begin(), chain.segs.end()));
  }

  return regions;
}

// ==========================================
// SegmentSelector implementation.
// ==========================================
namespace SegmentSelector {

std::vector<SegmentBool> select(const std::vector<SegmentBool> &segments,
                                const int selection[16]) {
  std::vector<SegmentBool> result;
  for (const auto &seg : segments) {
    int index =
        (seg.myFill.above.value_or(false) ? 8 : 0) +
        (seg.myFill.below.value_or(false) ? 4 : 0) +
        (seg.otherFill && seg.otherFill->above.value_or(false) ? 2 : 0) +
        (seg.otherFill && seg.otherFill->below.value_or(false) ? 1 : 0);

    int flags = selection[index];
    bool above = (flags & 1) != 0;
    bool below = (flags & 2) != 0;

    if ((!seg.closed && flags != 0) || (seg.closed && above != below)) {
      SegmentBool newSeg = seg;
      newSeg.myFill = {above, below};
      result.push_back(newSeg);
    }
  }
  return result;
}

std::vector<SegmentBool> union_(const std::vector<SegmentBool> &segments) {
  const int selection[16] = {4, 2, 1, 0, 2, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0};
  return select(segments, selection);
}

std::vector<SegmentBool> intersect(const std::vector<SegmentBool> &segments) {
  const int selection[16] = {0, 0, 0, 4, 0, 2, 0, 2, 0, 0, 1, 1, 4, 2, 1, 0};
  return select(segments, selection);
}

std::vector<SegmentBool> difference(const std::vector<SegmentBool> &segments) {
  const int selection[16] = {4, 0, 0, 0, 2, 0, 2, 0, 1, 1, 0, 0, 0, 1, 2, 0};
  return select(segments, selection);
}

std::vector<SegmentBool>
differenceRev(const std::vector<SegmentBool> &segments) {
  const int selection[16] = {4, 2, 1, 0, 0, 0, 1, 1, 0, 2, 0, 2, 0, 0, 0, 0};
  return select(segments, selection);
}

std::vector<SegmentBool> xor_(const std::vector<SegmentBool> &segments) {
  const int selection[16] = {4, 2, 1, 0, 2, 0, 0, 1, 1, 0, 0, 2, 0, 1, 2, 0};
  return select(segments, selection);
}

} // namespace SegmentSelector

// ==========================================
// Shape implementation.
// ==========================================
Shape::Shape(Shape &&) noexcept = default;
Shape &Shape::operator=(Shape &&) noexcept = default;
Shape::~Shape() = default;

Shape::Shape(const Geometry *geo,
             std::optional<std::vector<SegmentBool>> initSegments)
    : geo(geo) {
  if (initSegments.has_value()) {
    state = State::Seg;
    segmentsData = std::move(initSegments.value());
  } else {
    state = State::New;
    selfIntersecter = std::make_unique<Intersecter>(true, geo);
  }
}

Vec2 Shape::transformPoint(double x, double y) const {
  return {matrix.v[0] * x + matrix.v[2] * y + matrix.v[4],
          matrix.v[1] * x + matrix.v[3] * y + matrix.v[5]};
}

Shape &Shape::setTransform(double a, double b, double c, double d, double e,
                           double f) {
  if (state != State::New)
    throw std::runtime_error("PolyBool: Cannot change shape after using it");
  matrix = {a, b, c, d, e, f};
  return *this;
}

Shape &Shape::resetTransform() { return setTransform(1, 0, 0, 1, 0, 0); }
Vec6 Shape::getTransform() const {
  if (state != State::New)
    throw std::runtime_error("PolyBool: Cannot change shape after using it");
  return matrix;
}

Shape &Shape::transform(double a, double b, double c, double d, double e,
                        double f) {
  double a0 = matrix.v[0], b0 = matrix.v[1], c0 = matrix.v[2], d0 = matrix.v[3],
         e0 = matrix.v[4], f0 = matrix.v[5];
  matrix = {a0 * a + c0 * b, b0 * a + d0 * b,      a0 * c + c0 * d,
            b0 * c + d0 * d, a0 * e + c0 * f + e0, b0 * e + d0 * f + f0};
  return *this;
}

Shape &Shape::rotate(double angle) {
  double cos = std::cos(angle), sin = std::sin(angle);
  double a0 = matrix.v[0], b0 = matrix.v[1], c0 = matrix.v[2], d0 = matrix.v[3],
         e0 = matrix.v[4], f0 = matrix.v[5];
  matrix = {a0 * cos + c0 * sin,
            b0 * cos + d0 * sin,
            c0 * cos - a0 * sin,
            d0 * cos - b0 * sin,
            e0,
            f0};
  return *this;
}

Shape &Shape::rotateDeg(double angle) {
  double ang = std::fmod(std::fmod(angle, 360.0) + 360.0, 360.0);
  if (ang == 0.0)
    return *this;
  double rad = (PI * ang) / 180.0;
  return rotate(rad); // Keep the implementation compact by using the standard math path for all angles.
}

Shape &Shape::scale(double sx, double sy) {
  matrix.v[0] *= sx;
  matrix.v[1] *= sx;
  matrix.v[2] *= sy;
  matrix.v[3] *= sy;
  return *this;
}

Shape &Shape::translate(double tx, double ty) {
  matrix.v[4] += matrix.v[0] * tx + matrix.v[2] * ty;
  matrix.v[5] += matrix.v[1] * tx + matrix.v[3] * ty;
  return *this;
}

Shape &Shape::save() {
  if (state != State::New)
    throw std::runtime_error("PolyBool: Cannot change shape");
  saveStack.push_back(matrix);
  return *this;
}

Shape &Shape::restore() {
  if (state != State::New)
    throw std::runtime_error("PolyBool: Cannot change shape");
  if (!saveStack.empty()) {
    matrix = saveStack.back();
    saveStack.pop_back();
  }
  return *this;
}

Shape &Shape::beginPath() {
  if (state != State::New)
    throw std::runtime_error("PolyBool: Cannot change shape");
  selfIntersecter->beginPath();
  return endPath();
}

Shape &Shape::moveTo(double x, double y) {
  if (state != State::New)
    throw std::runtime_error("PolyBool: Cannot change shape");
  if (pathKind != PathKind::BeginPath)
    beginPath();
  Vec2 current = transformPoint(x, y);
  pathKind = PathKind::MoveTo;
  pathStart = current;
  pathCurrent = current;
  return *this;
}

Shape &Shape::lineTo(double x, double y) {
  if (state != State::New)
    throw std::runtime_error("PolyBool: Cannot change shape");
  if (pathKind != PathKind::MoveTo)
    throw std::runtime_error("PolyBool: Must call moveTo first");
  Vec2 current = transformPoint(x, y);
  selfIntersecter->addLine(pathCurrent, current);
  pathCurrent = current;
  return *this;
}

Shape &Shape::rect(double x, double y, double width, double height) {
  return moveTo(x, y)
      .lineTo(x + width, y)
      .lineTo(x + width, y + height)
      .lineTo(x, y + height)
      .closePath()
      .moveTo(x, y);
}

Shape &Shape::bezierCurveTo(double cp1x, double cp1y, double cp2x, double cp2y,
                            double x, double y) {
  if (state != State::New)
    throw std::runtime_error("PolyBool: Cannot change shape");
  if (pathKind != PathKind::MoveTo)
    throw std::runtime_error("PolyBool: Must call moveTo first");
  Vec2 current = transformPoint(x, y);
  selfIntersecter->addCurve(pathCurrent, transformPoint(cp1x, cp1y),
                            transformPoint(cp2x, cp2y), current);
  pathCurrent = current;
  return *this;
}

Shape &Shape::closePath() {
  if (state != State::New)
    throw std::runtime_error("PolyBool: Cannot change shape");
  if (pathKind == PathKind::MoveTo &&
      !geo->isEqualVec2(pathStart, pathCurrent)) {
    selfIntersecter->addLine(pathCurrent, pathStart);
    pathCurrent = pathStart;
  }
  selfIntersecter->closePath();
  return endPath();
}

Shape &Shape::endPath() {
  if (state != State::New)
    throw std::runtime_error("PolyBool: Cannot change shape");
  pathKind = PathKind::BeginPath;
  return *this;
}

const std::vector<SegmentBool> &Shape::selfIntersect() {
  if (state == State::New) {
    segmentsData = selfIntersecter->calculate();
    state = State::Seg;
    selfIntersecter.reset(); // Release the arena.
  }
  return segmentsData;
}

const std::vector<std::vector<Segment>> &Shape::segments() {
  if (state != State::Reg) {
    regionsData = SegmentChainer(selfIntersect(), geo);
    state = State::Reg;
  }
  return regionsData;
}

ShapeCombined Shape::combine(Shape &other) {
  Intersecter int_combine(false, geo);
  for (const auto &seg : selfIntersect()) {
    // Allocate through the Intersecter arena; it is cleaned up when int_combine is destroyed.
    SegmentBool *newSeg =
        int_combine.makeSegmentBool(seg.data, seg.myFill, seg.closed);
    int_combine.addSegment(newSeg, true);
  }
  for (const auto &seg : other.selfIntersect()) {
    SegmentBool *newSeg =
        int_combine.makeSegmentBool(seg.data, seg.myFill, seg.closed);
    int_combine.addSegment(newSeg, false);
  }
  return ShapeCombined(int_combine.calculate(), geo);
}

// ==========================================
// ShapeCombined implementation.
// ==========================================

Shape ShapeCombined::union_() const {
  return Shape(geo, SegmentSelector::union_(segments));
}
Shape ShapeCombined::intersect() const {
  return Shape(geo, SegmentSelector::intersect(segments));
}
Shape ShapeCombined::difference() const {
  return Shape(geo, SegmentSelector::difference(segments));
}
Shape ShapeCombined::differenceRev() const {
  return Shape(geo, SegmentSelector::differenceRev(segments));
}
Shape ShapeCombined::xor_() const {
  return Shape(geo, SegmentSelector::xor_(segments));
}

// ==========================================
// PolyBool API implementation.
// ==========================================

// Helper receiver for writing Shape output into Polygon regions.
class PolyPointReceiver : public IPolyBoolReceiver {
  std::vector<std::vector<PolyPoint>> *regions;

public:
  explicit PolyPointReceiver(std::vector<std::vector<PolyPoint>> *r)
      : regions(r) {}
  void beginPath() override {}
  void moveTo(double, double) override { regions->push_back({}); }
  void lineTo(double x, double y) override {
    regions->back().push_back(Vec2{x, y});
  }
  void bezierCurveTo(double c1x, double c1y, double c2x, double c2y, double x,
                     double y) override {
    regions->back().push_back(Vec6{c1x, c1y, c2x, c2y, x, y});
  }
  void closePath() override {}
};

Segments PolyBool::segments(const Polygon &poly) const {
  Shape s = shape();
  s.beginPath();
  for (const auto &region : poly.regions) {
    if (region.empty())
      continue;
    auto lastPoint = region.back();
    if (std::holds_alternative<Vec2>(lastPoint)) {
      s.moveTo(std::get<Vec2>(lastPoint).x, std::get<Vec2>(lastPoint).y);
    } else {
      s.moveTo(std::get<Vec6>(lastPoint).v[4], std::get<Vec6>(lastPoint).v[5]);
    }
    for (const auto &p : region) {
      if (std::holds_alternative<Vec2>(p)) {
        s.lineTo(std::get<Vec2>(p).x, std::get<Vec2>(p).y);
      } else {
        const auto &v6 = std::get<Vec6>(p);
        s.bezierCurveTo(v6.v[0], v6.v[1], v6.v[2], v6.v[3], v6.v[4], v6.v[5]);
      }
    }
    s.closePath();
  }
  return {std::move(s), poly.inverted};
}

CombinedSegments PolyBool::combine(Segments &segments1,
                                   Segments &segments2) const {
  return {segments1.shape.combine(segments2.shape), segments1.inverted,
          segments2.inverted};
}

Segments PolyBool::selectUnion(const CombinedSegments &combined) const {
  Shape s = combined.inverted1
                ? (combined.inverted2 ? combined.shape.intersect()
                                      : combined.shape.difference())
                : (combined.inverted2 ? combined.shape.differenceRev()
                                      : combined.shape.union_());
  return {std::move(s), combined.inverted1 || combined.inverted2};
}

Segments PolyBool::selectIntersect(const CombinedSegments &combined) const {
  Shape s = combined.inverted1
                ? (combined.inverted2 ? combined.shape.union_()
                                      : combined.shape.differenceRev())
                : (combined.inverted2 ? combined.shape.difference()
                                      : combined.shape.intersect());
  return {std::move(s), combined.inverted1 && combined.inverted2};
}

Segments PolyBool::selectDifference(const CombinedSegments &combined) const {
  Shape s = combined.inverted1
                ? (combined.inverted2 ? combined.shape.differenceRev()
                                      : combined.shape.union_())
                : (combined.inverted2 ? combined.shape.intersect()
                                      : combined.shape.difference());
  return {std::move(s), combined.inverted1 && !combined.inverted2};
}

Segments PolyBool::selectDifferenceRev(const CombinedSegments &combined) const {
  Shape s = combined.inverted1
                ? (combined.inverted2 ? combined.shape.difference()
                                      : combined.shape.intersect())
                : (combined.inverted2 ? combined.shape.union_()
                                      : combined.shape.differenceRev());
  return {std::move(s), !combined.inverted1 && combined.inverted2};
}

Segments PolyBool::selectXor(const CombinedSegments &combined) const {
  return {combined.shape.xor_(), combined.inverted1 != combined.inverted2};
}

Polygon PolyBool::polygon(Segments &segments) const {
  Polygon poly;
  poly.inverted = segments.inverted;
  PolyPointReceiver receiver(&poly.regions);
  segments.shape.output(receiver);
  return poly;
}

Polygon PolyBool::union_(const Polygon &poly1, const Polygon &poly2) const {
  auto seg1 = segments(poly1);
  auto seg2 = segments(poly2);
  auto comb = combine(seg1, seg2);
  auto seg3 = selectUnion(comb);
  return polygon(seg3);
}

Polygon PolyBool::intersect(const Polygon &poly1, const Polygon &poly2) const {
  auto seg1 = segments(poly1);
  auto seg2 = segments(poly2);
  auto comb = combine(seg1, seg2);
  auto seg3 = selectIntersect(comb);
  return polygon(seg3);
}

Polygon PolyBool::difference(const Polygon &poly1, const Polygon &poly2) const {
  auto seg1 = segments(poly1);
  auto seg2 = segments(poly2);
  auto comb = combine(seg1, seg2);
  auto seg3 = selectDifference(comb);
  return polygon(seg3);
}

Polygon PolyBool::differenceRev(const Polygon &poly1,
                                const Polygon &poly2) const {
  auto seg1 = segments(poly1);
  auto seg2 = segments(poly2);
  auto comb = combine(seg1, seg2);
  auto seg3 = selectDifferenceRev(comb);
  return polygon(seg3);
}

Polygon PolyBool::xor_(const Polygon &poly1, const Polygon &poly2) const {
  auto seg1 = segments(poly1);
  auto seg2 = segments(poly2);
  auto comb = combine(seg1, seg2);
  auto seg3 = selectXor(comb);
  return polygon(seg3);
}
} // namespace polybool
