#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <functional>
#include <optional>
#include <variant>
#include <vector>

namespace polybool {

template <typename T, size_t MaxCapacity> struct StackVector {
  T data[MaxCapacity];
  size_t count = 0;

  // Core std::vector-compatible interface.
  inline void push_back(const T &value) {
    if (count < MaxCapacity) {
      data[count++] = value;
    }
    // If MaxCapacity is exceeded, the value is ignored. Geometry code should
    // size these buffers so this path is not reached.
  }

  inline size_t size() const { return count; }
  inline bool empty() const { return count == 0; }
  inline void clear() { count = 0; }

  inline T &operator[](size_t i) { return data[i]; }
  inline const T &operator[](size_t i) const { return data[i]; }

  // Supports C++11 range-for loops: for (auto v : vec).
  inline T *begin() { return data; }
  inline T *end() { return data + count; }
  inline const T *begin() const { return data; }
  inline const T *end() const { return data + count; }
};

struct Vec2 {
  double x;
  double y;
};

struct Vec6 {
  double v[6];
};

struct BoundingBox {
  Vec2 min;
  Vec2 max;
};

// --- Standalone Math Helpers ---

inline double lerp(double a, double b, double t) { return a + (b - a) * t; }

inline Vec2 lerpVec2(const Vec2 &a, const Vec2 &b, double t) {
  return {lerp(a.x, b.x, t), lerp(a.y, b.y, t)};
}

inline bool boundingBoxesIntersect(const BoundingBox &bbox1,
                                   const BoundingBox &bbox2) {
  return !(bbox1.min.x > bbox2.max.x || bbox1.max.x < bbox2.min.x ||
           bbox1.min.y > bbox2.max.y || bbox1.max.y < bbox2.min.y);
}

// --- Geometry Interface Definitions ---
using Roots3 = StackVector<double, 3>;
using Roots2 = StackVector<double, 2>;

class Geometry {
public:
  virtual ~Geometry() = default;

  virtual double snap0(double v) const = 0;
  virtual double snap01(double v) const = 0;
  virtual bool isCollinear(const Vec2 &p1, const Vec2 &p2,
                           const Vec2 &p3) const = 0;
  virtual Roots3 solveCubic(double a, double b, double c, double d) const = 0;
  virtual bool isEqualVec2(const Vec2 &a, const Vec2 &b) const = 0;
  virtual int compareVec2(const Vec2 &a, const Vec2 &b) const = 0;
};

// Epsilon-based geometry implementation.
class GeometryEpsilon : public Geometry {
private:
  double epsilon;

  // Implemented internally in the cpp file.
  Roots3 solveCubicNormalized(double a, double b, double c) const;

public:
  explicit GeometryEpsilon(double eps = 1e-10) : epsilon(eps) {}

  double snap0(double v) const override {
    if (std::abs(v) < epsilon)
      return 0.0;
    return v;
  }

  double snap01(double v) const override {
    if (std::abs(v) < epsilon)
      return 0.0;
    if (std::abs(1.0 - v) < epsilon)
      return 1.0;
    return v;
  }

  bool isCollinear(const Vec2 &p1, const Vec2 &p2,
                   const Vec2 &p3) const override {
    double dx1 = p1.x - p2.x;
    double dy1 = p1.y - p2.y;
    double dx2 = p2.x - p3.x;
    double dy2 = p2.y - p3.y;
    return std::abs(dx1 * dy2 - dx2 * dy1) < epsilon;
  }

  bool isEqualVec2(const Vec2 &a, const Vec2 &b) const override {
    return std::abs(a.x - b.x) < epsilon && std::abs(a.y - b.y) < epsilon;
  }

  int compareVec2(const Vec2 &a, const Vec2 &b) const override {
    if (std::abs(b.x - a.x) < epsilon) {
      return std::abs(b.y - a.y) < epsilon ? 0 : (a.y < b.y ? -1 : 1);
    }
    return a.x < b.x ? -1 : 1;
  }

  // Complex root solving is implemented in the cpp file.
  Roots3 solveCubic(double a, double b, double c, double d) const override;
};

// --- Forward Declarations ---
class SegmentLine;
class SegmentCurve;
class IPolyBoolReceiver;

// --- Intersection Data Structures ---

struct TValuePair {
  double tA; // Progress on segment A (0.0 to 1.0).
  double tB; // Progress on segment B (0.0 to 1.0).
};

struct SegmentTValuePairs {
  StackVector<TValuePair, 9> tValuePairs; // Array of {seg1T, seg2T}.
};

struct SegmentTRangePairs {
  TValuePair tStart; // {seg1TStart, seg2TStart}
  TValuePair tEnd;   // {seg1TEnd, seg2TEnd}
};

// Use std::variant to model TypeScript SegmentTValuePairs |
// SegmentTRangePairs | null
using IntersectionResult =
    std::variant<std::monostate, SegmentTValuePairs, SegmentTRangePairs>;

// --- Builders ---

class SegmentTValuesBuilder {
public:
  StackVector<double, 8> tValues; // Replaces std::vector; stores up to 7 values.
  const Geometry *geo;

  explicit SegmentTValuesBuilder(const Geometry *geo) : geo(geo) {}

  template <typename Container>
  SegmentTValuesBuilder &addArray(const Container &ts) {
    for (double t : ts) {
      tValues.push_back(t); // StackVector::push_back guards capacity internally.
    }
    return *this;
  }

  SegmentTValuesBuilder &add(double t) {
    t = geo->snap01(t);
    // ignore values outside 0-1 range
    if (t < 0.0 || t > 1.0) {
      return *this;
    }
    for (double tv : tValues) {
      if (geo->snap0(t - tv) == 0.0) {
        // already have this location
        return *this;
      }
    }
    tValues.push_back(t);
    return *this;
  }

  StackVector<double, 8> list() {
    std::sort(tValues.begin(), tValues.end());
    return tValues;
  }
};

class SegmentTValuePairsBuilder {
public:
  StackVector<TValuePair, 9> tValuePairs;
  bool allowOutOfRange;
  const Geometry *geo;

  SegmentTValuePairsBuilder(bool allowOutOfRange, const Geometry *geo)
      : allowOutOfRange(allowOutOfRange), geo(geo) {}

  SegmentTValuePairsBuilder &add(double t1, double t2) {
    t1 = geo->snap01(t1);
    t2 = geo->snap01(t2);
    if (!allowOutOfRange && (t1 < 0.0 || t1 > 1.0 || t2 < 0.0 || t2 > 1.0))
      return *this;
    for (const auto &tv : tValuePairs) {
      if (geo->snap0(t1 - tv.tA) == 0.0 || geo->snap0(t2 - tv.tB) == 0.0)
        return *this;
    }
    tValuePairs.push_back({t1, t2});
    return *this;
  }

  StackVector<TValuePair, 9> list() {
    std::sort(tValuePairs.begin(), tValuePairs.end(),
              [](const auto &a, const auto &b) { return a.tA < b.tA; });
    return tValuePairs;
  }

  IntersectionResult done() {
    if (tValuePairs.empty())
      return std::monostate{};
    return SegmentTValuePairs{list()};
  }
};

// --- Segment Classes ---

using Segment = std::variant<SegmentLine, SegmentCurve>;

class SegmentLine {
public:
  Vec2 p0 = {0, 0};
  Vec2 p1 = {0, 0};
  const Geometry *geo = nullptr;

  SegmentLine() = default;

  SegmentLine(Vec2 p0, Vec2 p1, const Geometry *geo)
      : p0(p0), p1(p1), geo(geo) {}

  SegmentLine copy() const { return SegmentLine(p0, p1, geo); }
  bool isEqual(const SegmentLine &other) const {
    return geo->isEqualVec2(p0, other.p0) && geo->isEqualVec2(p1, other.p1);
  }
  Vec2 start() const { return p0; }
  Vec2 start2() const { return p1; }
  Vec2 end2() const { return p0; }
  Vec2 end() const { return p1; }

  void setStart(Vec2 p) { p0 = p; }
  void setEnd(Vec2 p) { p1 = p; }

  Vec2 point(double t) const {
    if (t == 0.0)
      return p0;
    if (t == 1.0)
      return p1;
    return {p0.x + (p1.x - p0.x) * t, p0.y + (p1.y - p0.y) * t};
  }

  template <typename OutputContainer>
  void split(const double *ts, size_t ts_count,
             OutputContainer &outResult) const {
    if (ts_count == 0) {
      outResult.push_back(*this);
      return;
    }
    Vec2 last = p0;
    for (size_t i = 0; i < ts_count; ++i) {
      Vec2 p = point(ts[i]);
      outResult.push_back(SegmentLine(last, p, geo));
      last = p;
    }
    outResult.push_back(SegmentLine(last, p1, geo));
  }

  template <typename OutputContainer>
  void splitSingle(double t, OutputContainer &outResult) const {
    split(&t, 1, outResult);
  }

  SegmentLine reverse() const { return SegmentLine(p1, p0, geo); }

  BoundingBox boundingBox() const {
    return {{std::min(p0.x, p1.x), std::min(p0.y, p1.y)},
            {std::max(p0.x, p1.x), std::max(p0.y, p1.y)}};
  }

  bool pointOn(const Vec2 &p) const { return geo->isCollinear(p, p0, p1); }

  template <typename TRecv> TRecv &draw(TRecv &ctx) const {
    static_assert(std::is_base_of_v<IPolyBoolReceiver, TRecv>,
                  "Receiver must inherit from IPolyBoolReceiver");

    ctx.moveTo(p0.x, p0.y);
    ctx.lineTo(p1.x, p1.y);
    return ctx;
  }
};

class SegmentCurve {
public:
  Vec2 p0 = {0, 0}, p1 = {0, 0}, p2 = {0, 0}, p3 = {0, 0}; // Default values.
  const Geometry *geo = nullptr;

  SegmentCurve() = default;

  SegmentCurve(Vec2 p0, Vec2 p1, Vec2 p2, Vec2 p3, const Geometry *geo)
      : p0(p0), p1(p1), p2(p2), p3(p3), geo(geo) {}

  SegmentCurve copy() const { return SegmentCurve(p0, p1, p2, p3, geo); }

  bool isEqual(const SegmentCurve &other) const {
    return geo->isEqualVec2(p0, other.p0) && geo->isEqualVec2(p1, other.p1) &&
           geo->isEqualVec2(p2, other.p2) && geo->isEqualVec2(p3, other.p3);
  }

  Vec2 start() const { return p0; }
  Vec2 start2() const { return p1; }
  Vec2 end2() const { return p2; }
  Vec2 end() const { return p3; }

  void setStart(Vec2 p) { p0 = p; }
  void setEnd(Vec2 p) { p3 = p; }

  Vec2 point(double t) const {
    if (t == 0.0)
      return p0;
    if (t == 1.0)
      return p3;

    double t1t = (1.0 - t) * (1.0 - t);
    double tt = t * t;
    double t0 = t1t * (1.0 - t);
    double t1 = 3.0 * t1t * t;
    double t2 = 3.0 * tt * (1.0 - t);
    double t3 = tt * t;

    return {p0.x * t0 + p1.x * t1 + p2.x * t2 + p3.x * t3,
            p0.y * t0 + p1.y * t1 + p2.y * t2 + p3.y * t3};
  }

  template <typename OutputContainer>
  void split(const double *ts, size_t ts_count,
             OutputContainer &outResult) const {
    if (ts_count == 0) {
      outResult.push_back(*this);
      return;
    }
    std::array<Vec2, 4> lastPts = {p0, p1, p2, p3};
    double lastT = 0.0;
    for (size_t i = 0; i < ts_count; ++i) {
      double t = ts[i];
      double localT = (t - lastT) / (1.0 - lastT);
      Vec2 p4 = lerpVec2(lastPts[0], lastPts[1], localT);
      Vec2 p5 = lerpVec2(lastPts[1], lastPts[2], localT);
      Vec2 p6 = lerpVec2(lastPts[2], lastPts[3], localT);
      Vec2 p7 = lerpVec2(p4, p5, localT);
      Vec2 p8 = lerpVec2(p5, p6, localT);
      Vec2 p9 = lerpVec2(p7, p8, localT);

      outResult.push_back(SegmentCurve(lastPts[0], p4, p7, p9, geo));
      lastPts = {p9, p8, p6, lastPts[3]};
      lastT = t;
    }
    outResult.push_back(
        SegmentCurve(lastPts[0], lastPts[1], lastPts[2], lastPts[3], geo));
  }

  template <typename OutputContainer>
  void splitSingle(double t, OutputContainer &outResult) const {
    split(&t, 1, outResult);
  }

  SegmentCurve reverse() const { return SegmentCurve(p3, p2, p1, p0, geo); }

  std::array<double, 4> getCubicCoefficients(int axis) const {
    double pp0 = (axis == 0) ? p0.x : p0.y;
    double pp1 = (axis == 0) ? p1.x : p1.y;
    double pp2 = (axis == 0) ? p2.x : p2.y;
    double pp3 = (axis == 0) ? p3.x : p3.y;
    return {pp3 - 3.0 * pp2 + 3.0 * pp1 - pp0,
            3.0 * pp2 - 6.0 * pp1 + 3.0 * pp0, 3.0 * pp1 - 3.0 * pp0, pp0};
  }

  // Defined in the cpp file.
  StackVector<double, 8> boundingTValues() const;
  StackVector<double, 8> inflectionTValues() const;
  BoundingBox boundingBox() const;

  std::optional<double> mapXtoT(double x, bool force = false) const;
  std::optional<double> mapXtoY(double x, bool force = false) const;
  bool pointOn(const Vec2 &p) const;

  std::optional<SegmentLine> toLine() const;

  template <typename TRecv> TRecv &draw(TRecv &ctx) const {
    static_assert(std::is_base_of_v<IPolyBoolReceiver, TRecv>,
                  "Receiver must inherit from IPolyBoolReceiver");

    ctx.moveTo(p0.x, p0.y);
    ctx.bezierCurveTo(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
    return ctx;
  }
};

// --- Intersection Function Declarations ---
double projectPointOntoSegmentLine(const Vec2 &p, const SegmentLine &seg);

IntersectionResult segmentLineIntersectSegmentLine(const SegmentLine &segA,
                                                   const SegmentLine &segB,
                                                   bool allowOutOfRange);
IntersectionResult segmentLineIntersectSegmentCurve(const SegmentLine &segA,
                                                    const SegmentCurve &segB,
                                                    bool allowOutOfRange,
                                                    bool invert);
IntersectionResult segmentCurveIntersectSegmentCurve(const SegmentCurve &segA,
                                                     const SegmentCurve &segB,
                                                     bool allowOutOfRange);
IntersectionResult segmentsIntersect(const Segment &segA, const Segment &segB,
                                     bool allowOutOfRange);

struct SegmentBoolFill {
  std::optional<bool> above = std::nullopt;
  std::optional<bool> below = std::nullopt;
};

template <typename T> struct ListBoolTransition {
  T before;
  T after;
  std::function<T(T)> insert;
};

struct SegmentBool {
  int id = -1;
  Segment data; // SegmentLine or SegmentCurve.
  SegmentBoolFill myFill;
  std::optional<SegmentBoolFill> otherFill = std::nullopt;
  bool closed = false;

  SegmentBool(const Segment &data, const SegmentBoolFill &fill = {},
              bool closed = false)
      : data(data), myFill(fill), closed(closed) {}
};

// C++ assignment performs deep copy here; copySegmentBool is not needed.

struct EventBool {
  bool isStart;
  Vec2 p;
  SegmentBool *seg; // Points to a SegmentBool owned by the arena.
  bool primary;
  EventBool *other = nullptr;
  EventBool *status = nullptr;

  EventBool(bool isStart, const Vec2 &p, SegmentBool *seg, bool primary)
      : isStart(isStart), p(p), seg(seg), primary(primary) {}
};

template <typename T> class ListBool {
public:
  // Array with a logical head pointer to keep removeHead inexpensive.
  std::vector<T> nodes;
  int headIndex = 0;

  void remove(T node) {
    auto it = std::find(nodes.begin(), nodes.end(), node);
    if (it != nodes.end())
      nodes.erase(it);
  }

  int getIndex(T node) const {
    auto it = std::find(nodes.begin(), nodes.end(), node);
    return it != nodes.end() ? std::distance(nodes.begin(), it) : -1;
  }

  bool isEmpty() const { return static_cast<size_t>(headIndex) >= nodes.size(); }
  T getHead() const { return nodes[headIndex]; }
  void removeHead() { headIndex++; }
  template <typename Func> void insertBefore(T node, Func &&check) {
    findTransition(node, std::forward<Func>(check)).insert(node);
  }

  template <typename Func>
  ListBoolTransition<T> findTransition(T node, Func &&check) {
    int i = headIndex;
    int high = nodes.size();
    while (i < high) {
      int mid = (i + high) >> 1;
      if (check(node) - check(nodes[mid]) > 0) {
        high = mid;
      } else {
        i = mid + 1;
      }
    }
    return {i <= headIndex ? nullptr : nodes[i - 1],
            i >= (int)nodes.size() ? nullptr : nodes[i], [this, i](T n) {
              this->nodes.insert(this->nodes.begin() + i, n);
              return n;
            }};
  }
};

// --- Intersecter Core Class Declaration ---
class Intersecter {
private:
  bool selfIntersection;
  const Geometry *geo;

  // Arena storage: appended elements keep stable pointers.
  std::deque<SegmentBool> segmentArena;
  std::deque<EventBool> eventArena;

  ListBool<EventBool *> events;
  ListBool<EventBool *> status;

  std::vector<SegmentBool *> currentPath;
  std::vector<Segment> splitBuffer;

public:

  // Factory methods: safely create objects inside the arena.
  SegmentBool *makeSegmentBool(const Segment &data,
                               const SegmentBoolFill &fill = {},
                               bool closed = false);
  EventBool *makeEventBool(bool isStart, const Vec2 &p, SegmentBool *seg,
                           bool primary);
  Intersecter(bool selfIntersection, const Geometry *geo)
      : selfIntersection(selfIntersection), geo(geo) {}

  int compareEvents(bool aStart, const Vec2 &a1, const Vec2 &a2,
                    const Segment &aSeg, bool bStart, const Vec2 &b1,
                    const Vec2 &b2, const Segment &bSeg);

  int compareSegments(const Segment &seg1, const Segment &seg2);

  void addEvent(EventBool *ev);
  EventBool *divideEvent(EventBool *ev, double t, const Vec2 &p);
  EventBool *addSegment(SegmentBool *seg, bool primary);

  void beginPath();
  void closePath();
  void addLine(const Vec2 &from, const Vec2 &to, bool primary = true);
  void addCurve(const Vec2 &from, const Vec2 &c1, const Vec2 &c2,
                const Vec2 &to, bool primary = true);

  ListBoolTransition<EventBool *> statusFindSurrounding(EventBool *ev);
  EventBool *checkIntersection(EventBool *ev1, EventBool *ev2);

  // Sweep-line main loop.
  std::vector<SegmentBool> calculate();
};

class IPolyBoolReceiver {
public:
  virtual ~IPolyBoolReceiver() = default;

  virtual void beginPath() = 0;
  virtual void moveTo(double x, double y) = 0;
  virtual void lineTo(double x, double y) = 0;
  virtual void bezierCurveTo(double cp1x, double cp1y, double cp2x, double cp2y,
                             double x, double y) = 0;
  virtual void closePath() = 0;
};

// --- 3. Merge and Chaining Function Declarations ---
std::optional<SegmentLine> joinLines(const SegmentLine &seg1,
                                     const SegmentLine &seg2,
                                     const Geometry *geo);
std::optional<SegmentCurve> joinCurves(const SegmentCurve &seg1,
                                       const SegmentCurve &seg2,
                                       const Geometry *geo);
std::optional<Segment> joinSegments(const Segment &seg1, const Segment &seg2,
                                    const Geometry *geo);

// Log parameter removed.
std::vector<std::vector<Segment>>
SegmentChainer(const std::vector<SegmentBool> &segments, const Geometry *geo);

// --- 4. Receiver Output Function Declaration ---
template <typename TRecv>
TRecv &segmentsToReceiver(const std::vector<std::vector<Segment>> &regions,
                          const Geometry *geo, TRecv &receiver,
                          const Vec6 &matrix) {
  static_assert(std::is_base_of_v<IPolyBoolReceiver, TRecv>,
                "Receiver must inherit from IPolyBoolReceiver");

  double a = matrix.v[0], b = matrix.v[1], c = matrix.v[2], d = matrix.v[3],
         e = matrix.v[4], f = matrix.v[5];

  receiver.beginPath();
  for (const auto &region : regions) {
    if (region.empty())
      continue;

    for (size_t i = 0; i < region.size(); ++i) {
      const auto &seg = region[i];

      // Use std::visit to handle the variant.
      std::visit(
          [&](auto &&arg) {
            using T = std::decay_t<decltype(arg)>;

            if (i == 0) {
              Vec2 p0 = arg.start();
              receiver.moveTo(a * p0.x + c * p0.y + e, b * p0.x + d * p0.y + f);
            }

            if constexpr (std::is_same_v<T, SegmentLine>) {
              Vec2 p1 = arg.p1;
              receiver.lineTo(a * p1.x + c * p1.y + e, b * p1.x + d * p1.y + f);
            } else if constexpr (std::is_same_v<T, SegmentCurve>) {
              Vec2 p1 = arg.p1;
              Vec2 p2 = arg.p2;
              Vec2 p3 = arg.p3;
              receiver.bezierCurveTo(
                  a * p1.x + c * p1.y + e, b * p1.x + d * p1.y + f,
                  a * p2.x + c * p2.y + e, b * p2.x + d * p2.y + f,
                  a * p3.x + c * p3.y + e, b * p3.x + d * p3.y + f);
            }
          },
          seg);
    }
    // Check whether the region is closed using std::visit.
    Vec2 startPt =
        std::visit([](auto &&arg) { return arg.start(); }, region.front());
    Vec2 endPt =
        std::visit([](auto &&arg) { return arg.end(); }, region.back());

    if (geo->isEqualVec2(startPt, endPt)) {
      receiver.closePath();
    }
  }
  return receiver;
}

// --- SegmentSelector Logic ---
namespace SegmentSelector {
std::vector<SegmentBool> select(const std::vector<SegmentBool> &segments,
                                const int selection[16]);
std::vector<SegmentBool> union_(const std::vector<SegmentBool> &segments);
std::vector<SegmentBool> intersect(const std::vector<SegmentBool> &segments);
std::vector<SegmentBool> difference(const std::vector<SegmentBool> &segments);
std::vector<SegmentBool>
differenceRev(const std::vector<SegmentBool> &segments);
std::vector<SegmentBool> xor_(const std::vector<SegmentBool> &segments);
} // namespace SegmentSelector

// --- High-Level API Data Structures ---
using PolyPoint = std::variant<Vec2, Vec6>;

struct Polygon {
  std::vector<std::vector<PolyPoint>> regions;
  bool inverted = false;
};

class ShapeCombined; // Forward declaration.

// --- Shape Path Builder ---
class Shape {
private:
  const Geometry *geo;

  enum class State { New, Seg, Reg };
  State state = State::New;

  // Use unique_ptr to avoid accidental Intersecter copies that would invalidate arena pointers.
  std::unique_ptr<Intersecter> selfIntersecter;
  std::vector<SegmentBool> segmentsData;
  std::vector<std::vector<Segment>> regionsData;

  enum class PathKind { BeginPath, MoveTo };
  PathKind pathKind = PathKind::BeginPath;
  Vec2 pathStart = {0, 0};
  Vec2 pathCurrent = {0, 0};

  std::vector<Vec6> saveStack;
  Vec6 matrix = {1, 0, 0, 1, 0, 0};

  Vec2 transformPoint(double x, double y) const;
  const std::vector<SegmentBool> &selfIntersect();

public:
  Shape(const Geometry *geo,
        std::optional<std::vector<SegmentBool>> initSegments = std::nullopt);

  // Movable, but not copyable.
  Shape(Shape &&) noexcept;
  Shape &operator=(Shape &&) noexcept;
  ~Shape();
  Shape(const Shape &) = delete;
  Shape &operator=(const Shape &) = delete;

  Shape &setTransform(double a, double b, double c, double d, double e,
                      double f);
  Shape &resetTransform();
  Vec6 getTransform() const;
  Shape &transform(double a, double b, double c, double d, double e, double f);
  Shape &rotate(double angle);
  Shape &rotateDeg(double angle);
  Shape &scale(double sx, double sy);
  Shape &translate(double tx, double ty);
  Shape &save();
  Shape &restore();

  Shape &beginPath();
  Shape &moveTo(double x, double y);
  Shape &lineTo(double x, double y);
  Shape &rect(double x, double y, double width, double height);
  Shape &bezierCurveTo(double cp1x, double cp1y, double cp2x, double cp2y,
                       double x, double y);
  Shape &closePath();
  Shape &endPath();

  const std::vector<std::vector<Segment>> &segments();

  template <typename TRecv>
  TRecv &output(TRecv &receiver, Vec6 outMatrix = {1, 0, 0, 1, 0, 0}) {
    return segmentsToReceiver(segments(), geo, receiver, outMatrix);
  }

  ShapeCombined combine(Shape &other);
};

// --- ShapeCombined Combiner ---
class ShapeCombined {
private:
  const Geometry *geo;
  std::vector<SegmentBool> segments;

public:
  ShapeCombined(const std::vector<SegmentBool> &segments, const Geometry *geo)
      : geo(geo), segments(segments) {}

  Shape union_() const;
  Shape intersect() const;
  Shape difference() const;
  Shape differenceRev() const;
  Shape xor_() const;
};

// --- PolyBool Top-Level Facade API ---
struct Segments {
  Shape shape;
  bool inverted;
};

struct CombinedSegments {
  ShapeCombined shape;
  bool inverted1;
  bool inverted2;
};

class PolyBool {
private:
  const Geometry *geo;

public:
  explicit PolyBool(const Geometry *geo) : geo(geo) {}

  Shape shape() const { return Shape(geo); }

  Segments segments(const Polygon &poly) const;
  CombinedSegments combine(Segments &segments1, Segments &segments2) const;

  Segments selectUnion(const CombinedSegments &combined) const;
  Segments selectIntersect(const CombinedSegments &combined) const;
  Segments selectDifference(const CombinedSegments &combined) const;
  Segments selectDifferenceRev(const CombinedSegments &combined) const;
  Segments selectXor(const CombinedSegments &combined) const;

  Polygon polygon(Segments &segments) const;

  // Convenience boolean operation API.
  Polygon union_(const Polygon &poly1, const Polygon &poly2) const;
  Polygon intersect(const Polygon &poly1, const Polygon &poly2) const;
  Polygon difference(const Polygon &poly1, const Polygon &poly2) const;
  Polygon differenceRev(const Polygon &poly1, const Polygon &poly2) const;
  Polygon xor_(const Polygon &poly1, const Polygon &poly2) const;
};

} // namespace polybool
