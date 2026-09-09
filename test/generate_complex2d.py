#!/usr/bin/env python3
"""
generate_complex2d.py - Generates complex 2D synthetic benchmark data for CEIF.

Geometric Components:
  1. Different scales: X spans [1000, 9000] (span ~7500), Y spans [0.2, 1.8] (span ~1.5) -> ~5000:1 aspect ratio.
  2. Main population with a hole: An annular donut structure centered at (5000, 1.0) with an interior void.
  3. Sharp spikes: Three narrow, needle-like rays extending outwards (North, North-East, South-West).
  4. Curved arc: A crescent / banana manifold hugging the perimeter from South to East.
  5. Separate spot: An isolated, compact Gaussian sub-population in the North-West quadrant.
  6. Ambient noise: A handful of sparse background points scattered across the bounding box.
"""

import math
import random
import sys


def generate_dataset(output_path='test/complex2d.csv', seed=42):
  random.seed(seed)
  points = []

  # 1. Main population with hole (donut: r in [2.2, 4.8])
  for _ in range(800):
    theta = random.uniform(0, 2 * math.pi)
    r = max(2.2, min(4.8, random.gauss(3.3, 0.45)))
    points.append((r * math.cos(theta), r * math.sin(theta)))

  # 2. Sharp spikes radiating outward
  # Spike 1: North (vertical needle)
  for _ in range(80):
    length = random.uniform(4.5, 9.0)
    u = random.gauss(0.0, 0.12)
    v = length
    points.append((u, v))

  # Spike 2: North-East (theta ~ 30 deg)
  for _ in range(80):
    length = random.uniform(4.5, 9.0)
    ang = math.radians(30) + random.gauss(0, 0.02)
    points.append((length * math.cos(ang), length * math.sin(ang)))

  # Spike 3: South-West (theta ~ 215 deg)
  for _ in range(80):
    length = random.uniform(4.5, 8.5)
    ang = math.radians(215) + random.gauss(0, 0.025)
    points.append((length * math.cos(ang), length * math.sin(ang)))

  # 3. Curved Arc (Crescent on South-East: r ~ 6.5, theta from -125 deg to 20 deg)
  for _ in range(320):
    theta = random.uniform(math.radians(-125), math.radians(20))
    r = random.gauss(6.5, 0.22)
    points.append((r * math.cos(theta), r * math.sin(theta)))

  # 4. Separate spot outside main population (North-West quadrant: isolated cluster)
  for _ in range(80):
    u = random.gauss(-7.0, 0.35)
    v = random.gauss(6.8, 0.35)
    points.append((u, v))

  # 5. Sparse ambient noise
  for _ in range(25):
    u = random.uniform(-9.5, 9.5)
    v = random.uniform(-9.5, 9.5)
    points.append((u, v))

  # Shuffle points to avoid ordering artifacts
  random.shuffle(points)

  # Scale coordinates: X in ~[1000, 9000], Y in ~[0.2, 1.8]
  with open(output_path, 'w') as f:
    for u, v in points:
      x = 5000.0 + u * 400.0
      y = 1.0 + v * 0.08
      f.write(f'{x:.4f},{y:.6f}\n')

  print(f'Wrote {len(points)} rows to {output_path}')


if __name__ == '__main__':
  path = sys.argv[1] if len(sys.argv) > 1 else 'test/complex2d.csv'
  generate_dataset(path)
