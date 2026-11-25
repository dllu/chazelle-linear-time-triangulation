"""Regenerate the README figures with python3 examples/gallery.py."""
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "docs" / "figures"


def snowflake():
    # Axial coordinates give exact 60-degree turns on a triangular lattice.
    ring = [(0, 0), (27, 0), (0, 27)]
    for _ in range(3):
        refined = []
        for (q, r), (s, t) in zip(ring, ring[1:] + ring[:1]):
            u, v = (s - q) // 3, (t - r) // 3
            refined.extend([(q, r), (q + u, r + v),
                            (q + 2*u + v, r + v - u), (q + 2*u, r + 2*v)])
        ring = refined
    return [(8 * (2*q + r), 14*r) for q, r in ring]


def hilbert_corridor():
    # Join four reflected / translated copies at each order.
    path = [(0, 0)]
    for width in [1, 2, 4, 8]:
        path = ([(y, x) for x, y in path]
                + [(x, y + width) for x, y in path]
                + [(x + width, y + width) for x, y in path]
                + [(2*width - 1 - y, width - 1 - x) for x, y in path])

    # Offset the two walls at each turn. A separation of 32 between
    # centerlines leaves a gap of 12 between neighboring passages.
    normals = [(y - v, u - x) for (x, y), (u, v) in zip(path, path[1:])]
    left, right = [], []
    for (x, y), before, after in zip(path, [(0, 0)] + normals, normals + [(0, 0)]):
        if before == after:
            continue
        dx = 10 * (before[0] + after[0])
        dy = 10 * (before[1] + after[1])
        left.append((32*x + dx, 32*y + dy))
        right.append((32*x - dx, 32*y - dy))
    return left + right[::-1]


def render(name, title, ring, fill, ink):
    polygon = ''.join(f'{x} {y}\n' for x, y in ring)
    run = subprocess.run(['cargo', 'run', '--quiet', '--release', '--bin', 'chazelle'],
                         cwd=ROOT, input=polygon, text=True, stdout=subprocess.PIPE, check=True)
    triangles = [tuple(map(int, line.split())) for line in run.stdout.splitlines()]
    assert len(triangles) == len(ring) - 2

    xs, ys = zip(*ring)
    side = max(max(xs) - min(xs), max(ys) - min(ys)) + 30
    left = (min(xs) + max(xs) - side) // 2
    top = (-min(ys) - max(ys) - side) // 2
    points = [f'{x},{-y}' for x, y in ring]  # SVG y increases down the page.
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="420" height="420" viewBox="{left} {top} {side} {side}" role="img">',
           f'<title>{title}</title>',
           f'<desc>{len(ring)} boundary vertices and {len(triangles)} triangles.</desc>',
           f'<rect x="{left}" y="{top}" width="{side}" height="{side}" rx="10" fill="#faf9f6"/>',
           f'<g fill="{fill}" stroke="{ink}" stroke-width="0.8" stroke-opacity="0.75" stroke-linejoin="round">']
    for a, b, c in triangles:
        (x, y), (u, v), (s, t) = ring[a], ring[b], ring[c]
        area2 = (u - x)*(t - y) - (v - y)*(s - x)
        assert area2 > 0
        opacity = min(0.8, 0.28 + 0.08 * max(0, area2.bit_length() - 8))
        svg.append(f'<polygon points="{points[a]} {points[b]} {points[c]}" fill-opacity="{opacity:.2f}"/>')
    boundary = ' '.join(points)
    svg.extend(['</g>', f'<polygon points="{boundary}" fill="none" stroke="{ink}" '
                'stroke-width="1.9" stroke-linejoin="round"/>', '</svg>'])
    (FIGURES / f'{name}.svg').write_text('\n'.join(svg) + '\n', encoding='utf-8')


if __name__ == '__main__':
    render('koch-snowflake', 'Koch snowflake triangulation', snowflake(), '#4ea59e', '#204b50')
    render('hilbert-corridor', 'Hilbert corridor triangulation', hilbert_corridor(), '#db9a5b', '#714227')
