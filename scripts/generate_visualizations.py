#!/usr/bin/env python3

from __future__ import annotations

import html
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "visualizations"

WIDTH = 430
HEIGHT = 560
PANEL_X = 32
PANEL_W = WIDTH - 64
TOP_Y = 24
BOTTOM_Y = 328
PANEL_H = 182

BLUE = "#1f38ff"
RED = "#ff2a2a"
GREEN = "#12b81f"
TEXT = "#202124"


def poly(points):
    return {"regions": [points], "inverted": False}


def multi(regions, inverted=False):
    return {"regions": regions, "inverted": inverted}


def empty(inverted=False):
    return {"regions": [], "inverted": inverted}


def curve_path(commands, inverted=False):
    return {"paths": [commands], "regions": [], "inverted": inverted}


def shape_path(regions, inverted=False):
    return {"regions": regions, "inverted": inverted}


def inverted(shape):
    copy = dict(shape)
    copy["inverted"] = True
    return copy


triangle1 = poly([(0, 0), (5, 10), (10, 0)])
triangle2 = poly([(5, 0), (10, 10), (15, 0)])
triangle_intersection = poly([(10, 0), (5, 0), (7.5, 5)])
triangle_union = poly([(10, 10), (7.5, 5), (5, 10), (0, 0), (15, 0)])

box1 = poly([(0, 0), (5, 0), (5, -5), (0, -5)])
curve1 = curve_path([("M", 0, 0), ("C", 0, -5, 10, -5, 10, 0), ("Z",)])
curve_union = curve_path(
    [
        ("M", 0, 0),
        ("L", 10, 0),
        ("C", 10, -2.5, 7.5, -3.75, 5, -3.75),
        ("L", 5, -5),
        ("L", 0, -5),
        ("Z",),
    ]
)

square_a = poly([(0, 0), (10, 0), (10, 10), (0, 10)])
square_overlap_b = poly([(5, 5), (15, 5), (15, 15), (5, 15)])
square_disjoint_b = poly([(20, 20), (30, 20), (30, 30), (20, 30)])
square_edge_touch_b = poly([(10, 0), (20, 0), (20, 10), (10, 10)])
square_point_touch_b = poly([(10, 10), (20, 10), (20, 20), (10, 20)])
square_inside_b = poly([(2, 2), (8, 2), (8, 8), (2, 8)])
horizontal_cutter = poly([(5, 0), (25, 0), (25, 10), (5, 10)])
two_squares = multi(
    [
        [(0, 0), (10, 0), (10, 10), (0, 10)],
        [(20, 0), (30, 0), (30, 10), (20, 10)],
    ]
)

overlap_union = poly(
    [(0, 0), (10, 0), (10, 5), (15, 5), (15, 15), (5, 15), (5, 10), (0, 10)]
)
overlap_difference = poly([(0, 0), (10, 0), (10, 5), (5, 5), (5, 10), (0, 10)])
overlap_difference_rev = poly(
    [(5, 10), (10, 10), (10, 5), (15, 5), (15, 15), (5, 15)]
)
overlap_xor = multi(
    [
        [(0, 0), (10, 0), (10, 5), (5, 5), (5, 10), (0, 10)],
        [(5, 10), (10, 10), (10, 5), (15, 5), (15, 15), (5, 15)],
    ]
)
contained_difference = multi(
    [
        [(0, 0), (10, 0), (10, 10), (0, 10)],
        [(2, 2), (8, 2), (8, 8), (2, 8)],
    ]
)
disjoint_union = multi(
    [
        [(0, 0), (10, 0), (10, 10), (0, 10)],
        [(20, 20), (30, 20), (30, 30), (20, 30)],
    ]
)
edge_union = poly([(0, 0), (20, 0), (20, 10), (0, 10)])
point_touch_union = multi(
    [
        [(0, 0), (10, 0), (10, 10), (0, 10)],
        [(10, 10), (20, 10), (20, 20), (10, 20)],
    ]
)
multi_intersection = multi(
    [
        [(5, 0), (10, 0), (10, 10), (5, 10)],
        [(20, 0), (25, 0), (25, 10), (20, 10)],
    ]
)
multi_difference = multi(
    [
        [(0, 0), (5, 0), (5, 10), (0, 10)],
        [(25, 0), (30, 0), (30, 10), (25, 10)],
    ]
)
inverted_union = multi([[(0, 0), (10, 0), (10, 5), (5, 5), (5, 10), (0, 10)]], True)
inverted_difference = multi(
    [[(0, 0), (10, 0), (10, 5), (15, 5), (15, 15), (5, 15), (5, 10), (0, 10)]],
    True,
)
inverted_xor = multi(
    [
        [(0, 0), (10, 0), (10, 5), (5, 5), (5, 10), (0, 10)],
        [(5, 10), (10, 10), (10, 5), (15, 5), (15, 15), (5, 15)],
    ],
    True,
)

shape1 = multi(
    [
        [(50, 50), (150, 150), (190, 50)],
        [(130, 50), (290, 150), (290, 50)],
    ]
)
shape2 = multi(
    [
        [(110, 20), (110, 110), (20, 20)],
        [(130, 170), (130, 20), (260, 20), (260, 170)],
    ]
)
shape_intersection = multi(
    [
        [(110, 110), (50, 50), (110, 50)],
        [(150, 150), (178, 80), (130, 50), (130, 130)],
        [(260, 131.25), (178, 80), (190, 50), (260, 50)],
    ]
)
transform_input = poly([(50, 50), (-10, 50), (10, 10)])
transform_output = poly([(250, 300), (70, 300), (130, 220)])


CASES = [
    ("triangles / intersect", "Intersect", triangle1, triangle2, triangle_intersection),
    ("triangles / union", "Union", triangle1, triangle2, triangle_union),
    ("curve / union with box", "Union", box1, curve1, curve_union),
    ("rect overlap / intersect", "Intersect", square_a, square_overlap_b, poly([(5, 5), (10, 5), (10, 10), (5, 10)])),
    ("rect overlap / union", "Union", square_a, square_overlap_b, overlap_union),
    ("rect overlap / difference", "Difference", square_a, square_overlap_b, overlap_difference),
    ("rect overlap / differenceRev", "DifferenceRev", square_a, square_overlap_b, overlap_difference_rev),
    ("rect overlap / xor", "Xor", square_a, square_overlap_b, overlap_xor),
    ("disjoint rects / intersect", "Intersect", square_a, square_disjoint_b, empty()),
    ("disjoint rects / union", "Union", square_a, square_disjoint_b, disjoint_union),
    ("disjoint rects / difference", "Difference", square_a, square_disjoint_b, square_a),
    ("identical rects / union", "Union", square_a, square_a, square_a),
    ("identical rects / intersect", "Intersect", square_a, square_a, square_a),
    ("identical rects / difference", "Difference", square_a, square_a, empty()),
    ("identical rects / xor", "Xor", square_a, square_a, empty()),
    ("shared edge rects / union", "Union", square_a, square_edge_touch_b, edge_union),
    ("shared edge rects / intersect", "Intersect", square_a, square_edge_touch_b, empty()),
    ("shared edge rects / xor", "Xor", square_a, square_edge_touch_b, edge_union),
    ("point-touch rects / union", "Union", square_a, square_point_touch_b, point_touch_union),
    ("point-touch rects / intersect", "Intersect", square_a, square_point_touch_b, empty()),
    ("contained rect / union", "Union", square_a, square_inside_b, square_a),
    ("contained rect / intersect", "Intersect", square_a, square_inside_b, square_inside_b),
    ("contained rect / difference", "Difference", square_a, square_inside_b, contained_difference),
    ("contained rect / xor", "Xor", square_a, square_inside_b, contained_difference),
    ("contained rect / differenceRev", "DifferenceRev", square_a, square_inside_b, empty()),
    ("multi-region / intersect", "Intersect", two_squares, horizontal_cutter, multi_intersection),
    ("multi-region / difference", "Difference", two_squares, horizontal_cutter, multi_difference),
    ("inverted first polygon / union", "Union", inverted(square_a), square_overlap_b, inverted_union),
    ("inverted first polygon / intersect", "Intersect", inverted(square_a), square_overlap_b, overlap_difference_rev),
    ("inverted first polygon / difference", "Difference", inverted(square_a), square_overlap_b, inverted_difference),
    ("inverted first polygon / xor", "Xor", inverted(square_a), square_overlap_b, inverted_xor),
    ("shape output / intersect example", "Intersect", shape1, shape2, shape_intersection),
    ("shape output / transform", "Transform", transform_input, None, transform_output),
]


def slugify(name):
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def all_points(shape):
    if shape is None:
        return []

    points = []
    for region in shape.get("regions", []):
        points.extend(region)

    for path in shape.get("paths", []):
        for cmd in path:
            if cmd[0] == "M" or cmd[0] == "L":
                points.append((cmd[1], cmd[2]))
            elif cmd[0] == "C":
                points.extend([(cmd[1], cmd[2]), (cmd[3], cmd[4]), (cmd[5], cmd[6])])
    return points


def bounds_for(shapes):
    points = []
    for shape in shapes:
        points.extend(all_points(shape))

    if not points:
        return (0.0, 0.0, 1.0, 1.0)

    min_x = min(p[0] for p in points)
    max_x = max(p[0] for p in points)
    min_y = min(p[1] for p in points)
    max_y = max(p[1] for p in points)
    if math.isclose(min_x, max_x):
        min_x -= 1
        max_x += 1
    if math.isclose(min_y, max_y):
        min_y -= 1
        max_y += 1

    pad_x = (max_x - min_x) * 0.12
    pad_y = (max_y - min_y) * 0.12
    return (min_x - pad_x, min_y - pad_y, max_x + pad_x, max_y + pad_y)


def transformer(bounds, y0, panel_h):
    min_x, min_y, max_x, max_y = bounds
    scale = min(PANEL_W / (max_x - min_x), panel_h / (max_y - min_y))
    used_w = (max_x - min_x) * scale
    used_h = (max_y - min_y) * scale
    offset_x = PANEL_X + (PANEL_W - used_w) / 2
    offset_y = y0 + (panel_h - used_h) / 2

    def tx(point):
        x, y = point
        return (offset_x + (x - min_x) * scale, offset_y + (y - min_y) * scale)

    return tx


def format_num(value):
    if abs(value - round(value)) < 1e-9:
        return str(int(round(value)))
    return f"{value:.3f}".rstrip("0").rstrip(".")


def region_path(region, tx):
    if not region:
        return ""
    x, y = tx(region[0])
    chunks = [f"M {format_num(x)} {format_num(y)}"]
    for point in region[1:]:
        x, y = tx(point)
        chunks.append(f"L {format_num(x)} {format_num(y)}")
    chunks.append("Z")
    return " ".join(chunks)


def command_path(commands, tx):
    chunks = []
    for cmd in commands:
        if cmd[0] == "M":
            x, y = tx((cmd[1], cmd[2]))
            chunks.append(f"M {format_num(x)} {format_num(y)}")
        elif cmd[0] == "L":
            x, y = tx((cmd[1], cmd[2]))
            chunks.append(f"L {format_num(x)} {format_num(y)}")
        elif cmd[0] == "C":
            x1, y1 = tx((cmd[1], cmd[2]))
            x2, y2 = tx((cmd[3], cmd[4]))
            x3, y3 = tx((cmd[5], cmd[6]))
            chunks.append(
                "C "
                f"{format_num(x1)} {format_num(y1)} "
                f"{format_num(x2)} {format_num(y2)} "
                f"{format_num(x3)} {format_num(y3)}"
            )
        elif cmd[0] == "Z":
            chunks.append("Z")
    return " ".join(chunks)


def shape_path_data(shape, tx, frame=None):
    pieces = []
    if shape is None:
        return ""

    if shape.get("inverted") and frame:
        x0, y0, x1, y1 = frame
        pieces.append(f"M {x0} {y0} L {x1} {y0} L {x1} {y1} L {x0} {y1} Z")

    for region in shape.get("regions", []):
        pieces.append(region_path(region, tx))
    for path in shape.get("paths", []):
        pieces.append(command_path(path, tx))
    return " ".join(piece for piece in pieces if piece)


def draw_vertices(shape, tx, color):
    if shape is None:
        return ""
    circles = []
    for point in all_points(shape):
        x, y = tx(point)
        circles.append(
            f'<circle cx="{format_num(x)}" cy="{format_num(y)}" r="3" fill="{color}"/>'
        )
    return "\n".join(circles)


def draw_shape(shape, tx, stroke, fill, opacity, frame):
    data = shape_path_data(shape, tx, frame)
    if not data:
        return ""
    fill_rule = "evenodd" if shape and shape.get("inverted") else "evenodd"
    return (
        f'<path d="{data}" fill="{fill}" fill-opacity="{opacity}" '
        f'stroke="{stroke}" stroke-width="1.5" fill-rule="{fill_rule}"/>'
    )


def arrow():
    cx = WIDTH / 2
    return (
        f'<line x1="{cx}" y1="257" x2="{cx}" y2="292" stroke="#000" stroke-width="13"/>'
        f'<polygon points="{cx - 18},286 {cx + 18},286 {cx},312" fill="#000"/>'
    )


def render_case(index, name, op, shape_a, shape_b, result):
    bounds = bounds_for([shape_a, shape_b, result])
    top_tx = transformer(bounds, TOP_Y, PANEL_H)
    bottom_tx = transformer(bounds, BOTTOM_Y, PANEL_H)
    frame_top = (PANEL_X, TOP_Y, PANEL_X + PANEL_W, TOP_Y + PANEL_H)
    frame_bottom = (PANEL_X, BOTTOM_Y, PANEL_X + PANEL_W, BOTTOM_Y + PANEL_H)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        '<rect width="100%" height="100%" fill="white"/>',
        draw_shape(shape_a, top_tx, BLUE, BLUE, "0.23", frame_top),
    ]

    if shape_b is not None:
        parts.append(draw_shape(shape_b, top_tx, RED, RED, "0.23", frame_top))

    parts.extend(
        [
            draw_vertices(shape_a, top_tx, BLUE),
            draw_vertices(shape_b, top_tx, RED),
            f'<text x="{WIDTH / 2}" y="240" text-anchor="middle" font-family="Arial, sans-serif" font-size="20" fill="{TEXT}">{html.escape(op)}</text>',
            arrow(),
        ]
    )

    if not all_points(result):
        parts.append(
            f'<text x="{WIDTH / 2}" y="{BOTTOM_Y + PANEL_H / 2 + 7}" text-anchor="middle" '
            f'font-family="Arial, sans-serif" font-size="22" fill="{GREEN}">Empty</text>'
        )
    else:
        parts.append(draw_shape(result, bottom_tx, GREEN, "#31ff31", "0.80", frame_bottom))
        parts.append(draw_vertices(result, bottom_tx, "#087d12"))

    parts.extend(
        [
            f'<text x="16" y="548" font-family="Arial, sans-serif" font-size="12" fill="#6b7280">{index:02d}. {html.escape(name)}</text>',
            "</svg>",
        ]
    )
    return "\n".join(part for part in parts if part)


def write_index(items):
    cards = []
    for index, name, filename in items:
        cards.append(
            f"""
      <a class="card" href="{filename}">
        <img src="{filename}" alt="{html.escape(name)}">
        <div><span>{index:02d}</span>{html.escape(name)}</div>
      </a>"""
        )

    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>polybool_cpp Test Visualizations</title>
  <style>
    body {{
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: #202124;
      background: #f6f7f9;
    }}
    header {{
      padding: 28px 32px 18px;
      background: #fff;
      border-bottom: 1px solid #e5e7eb;
    }}
    h1 {{
      margin: 0 0 8px;
      font-size: 24px;
      font-weight: 650;
    }}
    p {{
      margin: 0;
      color: #5f6368;
      line-height: 1.45;
    }}
    main {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(250px, 1fr));
      gap: 18px;
      padding: 24px 32px 40px;
    }}
    .card {{
      display: block;
      overflow: hidden;
      color: inherit;
      text-decoration: none;
      background: #fff;
      border: 1px solid #e5e7eb;
      border-radius: 8px;
      box-shadow: 0 1px 2px rgba(0, 0, 0, 0.04);
    }}
    .card img {{
      display: block;
      width: 100%;
      background: #fff;
      border-bottom: 1px solid #edf0f2;
    }}
    .card div {{
      padding: 10px 12px 12px;
      font-size: 14px;
    }}
    .card span {{
      display: inline-block;
      min-width: 28px;
      color: #6b7280;
    }}
  </style>
</head>
<body>
  <header>
    <h1>polybool_cpp Test Visualizations</h1>
    <p>Blue and red shapes are inputs. Green shapes are the expected results.</p>
  </header>
  <main>
{''.join(cards)}
  </main>
</body>
</html>
"""
    (OUT_DIR / "index.html").write_text(content)


def main():
    OUT_DIR.mkdir(exist_ok=True)
    items = []
    for index, (name, op, shape_a, shape_b, result) in enumerate(CASES, start=1):
        filename = f"{index:02d}_{slugify(name)}.svg"
        (OUT_DIR / filename).write_text(render_case(index, name, op, shape_a, shape_b, result))
        items.append((index, name, filename))
    write_index(items)
    print(f"Wrote {len(items)} SVG files and index.html to {OUT_DIR}")


if __name__ == "__main__":
    main()
