#!/usr/bin/env python

# Copyright Bohan Cao (2110313@mail.nankai.edu.cn) 2023
# Draw a svg figure of 3 dimensional triangular lattice with (m, n)

import drawsvg as dw
import numpy as np

m = 13
n = 7

fig = dw.Drawing(width=512, height=512, origin=(-40, -250))

A = np.arctan(np.sqrt(3)/2 * m / (n - m/2)) if 2*n != m else np.pi / 2
B = np.arctan(np.sqrt(3)/2 * n / (m - n/2)) if 2*m != n else np.pi / 2

k = 35

def pos(i, j):
    return -(i-n)*np.cos(A) + j*np.cos(B), (i-n)*np.sin(A) + j*np.sin(B)

for i in range(n+1):
    x, y = pos(i, 0)
    x2, y2 = pos(i, m)
    fig.append(dw.Line(k*x, -k*y, k*x2, -k*y2, stroke='black'))

for i in range(m+1):
    x, y = pos(0, i)
    x2, y2 = pos(n, i)
    fig.append(dw.Line(k*x, -k*y, k*x2, -k*y2, stroke='black'))

for i in range(m+n+1):
    x, y = pos(i, 0) if i<n else pos(n, i-n)
    x2, y2 = pos(0, i) if i<m else pos(i-m, m)
    fig.append(dw.Line(k*x, -k*y, k*x2, -k*y2, stroke='black'))

fig.save_svg('grid_3.svg')
