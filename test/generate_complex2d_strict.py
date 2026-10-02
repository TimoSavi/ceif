#!/usr/bin/env python3
"""
generate_complex2d_strict.py - Generates complex 2D synthetic benchmark data
with strict, razor-sharp geometric boundaries and ZERO outliers.

Geometric Components (Supervised / Clean Inliers):
  1. Main population with strict hole: Uniform annular donut with r in [2.2, 4.5] (hard boundaries).
  2. Sharp rectangular spikes:
     - North: u in [-0.20, 0.20], v in [4.5, 9.0]
     - North-East (30 deg): length in [4.5, 9.0], width in [-0.20, 0.20]
     - South-West (215 deg): length in [4.5, 8.5], width in [-0.20, 0.20]
  3. Curved arc: Circular sector band r in [6.2, 6.8], theta in [-125 deg, 20 deg]
  4. Separate spot: Hard circular disk centered at (-7.0, 6.8) with radius R = 0.8
  5. Ambient noise: 0 points (100% pure inliers, zero contaminants).
"""

import math
import random
import sys


def generate_dataset(output_path='test/complex2d_strict.csv', seed=42):
  random.seed(seed)
  points = []

  # 1. Main population with hole (strict donut: r in [2.2, 4.5], area-uniform)
  r_min_sq = 2.2**2
  r_max_sq = 4.5**2
  for _ in range(800):
    theta = random.uniform(0, 2 * math.pi)
    r = math.sqrt(random.uniform(r_min_sq, r_max_sq))
    points.append((r * math.cos(theta), r * math.sin(theta)))

  # 2. Sharp spikes radiating outward (strict rectangular strips)
  # Spike 1: North (vertical strip)
  for _ in range(80):
    u = random.uniform(-0.20, 0.20)
    v = random.uniform(4.5, 9.0)
    points.append((u, v))

  # Spike 2: North-East (theta = 30 deg, hard rectangular strip)
  ang2 = math.radians(30)
  cos2, sin2 = math.cos(ang2), math.sin(ang2)
  for _ in range(80):
    length = random.uniform(4.5, 9.0)
    width = random.uniform(-0.20, 0.20)
    u = length * cos2 - width * sin2
    v = length * sin2 + width * cos2
    points.append((u, v))

  # Spike 3: South-West (theta = 215 deg, hard rectangular strip)
  ang3 = math.radians(215)
  cos3, sin3 = math.cos(ang3), math.sin(ang3)
  for _ in range(80):
    length = random.uniform(4.5, 8.5)
    width = random.uniform(-0.20, 0.20)
    u = length * cos3 - width * sin3
    v = length * sin3 + width * cos3
    points.append((u, v))

  # 3. Curved Arc (Crescent on South-East: r in [6.2, 6.8], theta in [-125 deg, 20 deg])
  arc_r_min_sq = 6.2**2
  arc_r_max_sq = 6.8**2
  theta_start = math.radians(-125)
  theta_end = math.radians(20)
  for _ in range(320):
    theta = random.uniform(theta_start, theta_end)
    r = math.sqrt(random.uniform(arc_r_min_sq, arc_r_max_sq))
    points.append((r * math.cos(theta), r * math.sin(theta)))

  # 4. Separate spot outside main population (Hard circular disk at (-7.0, 6.8), R = 0.8)
  for _ in range(80):
    r = 0.8 * math.sqrt(random.uniform(0.0, 1.0))
    phi = random.uniform(0, 2 * math.pi)
    u = -7.0 + r * math.cos(phi)
    v = 6.8 + r * math.sin(phi)
    points.append((u, v))

  # 5. Zero ambient noise!

  # Shuffle points to avoid ordering artifacts
  random.shuffle(points)

  # Scale coordinates: X in ~[1000, 9000], Y in ~[0.2, 1.8]
  with open(output_path, 'w') as f:
    for u, v in points:
      x = 5000.0 + u * 400.0
      y = 1.0 + v * 0.08
      f.write(f'{x:.4f},{y:.6f}\n')

  print(f'Wrote {len(points)} pure inlier rows with strict boundaries to {output_path}')


if __name__ == '__main__':
  path = sys.argv[1] if len(sys.argv) > 1 else 'test/complex2d_strict.csv'
  generate_dataset(path)
