#!/usr/bin/env python

# Force-Directed Graph Drawing

import Tkinter
import random
import math

import numpy as np

x_fixed = np.array([
    [ 0.747264, 90.93342 ],
    [ 0.806094, 91.614726],
    [ 0.741519, 88.184657],
    [ 0.723471, 88.817481],
    [ 0.777587, 90.1086  ],
    [ 0.725995, 87.924487],
    [ 0.834896, 91.501665],
    [ 0.722248, 88.508865],
    [ 0.767785, 88.731781]
])

x = np.array([
    [ 0.72241866, 91.75850967],
    [ 0.79773418, 92.75932702],
    [ 0.75105485, 87.05160468],
    [ 0.68909309, 89.00189335],
    [ 0.81716576, 88.80893972],
    [ 0.70543842, 86.97344667],
    [ 0.86375709, 90.83958238],
    [ 0.68778784, 88.35101251],
    [ 0.79267796, 87.90889185]
])

def norm(pts, mean=None, scaling=None):
    if mean is None:
        mean = pts.mean(axis=0)
    pts -= mean

    if scaling is None:
        scaling = 4*np.max(np.abs(x), axis=0)
    pts = pts/scaling + 0.5

    return pts, mean, scaling

x, pts_mean, pts_scaling = norm(x)
x_fixed, _, _ = norm(x_fixed, mean=pts_mean, scaling=pts_scaling)

v = np.zeros_like(x)


# mass
alpha = 1.0
beta = .0001
k = 0.2
#damping
eta = .99
delta_t = .01 
m = len(x_fixed)

dij = 0.1

UPDATE_FREQ = 1


root = Tkinter.Tk()
canvas = Tkinter.Canvas(root, width=500, height=500, background="yellow")
canvas.pack()

ids = []


def move_oval(id, xi):
    newx = int(xi[0] * 500)
    newy = int(xi[1] * 500)
    canvas.coords(id, newx - 5, newy - 5, newx + 5, newy + 5)

for i in xrange(m):
    xi = x_fixed[i]
    id = canvas.create_oval(245, 245, 255, 255, fill="red")
    move_oval(id, xi)

    xi = x[i]
    id = canvas.create_oval(245, 245, 255, 255, fill="blue")
    ids.append(id)
    move_oval(id, xi)

lids = []


def move_line(id, xi, xj):
    canvas.coords(id,
                  int(xi[0] * 500),
                  int(xi[1] * 500),
                  int(xj[0] * 500),
                  int(xj[1] * 500))

for i in xrange(m):
    id = canvas.create_line(0, 0, 0, 0)
    lids.append(id)
    move_line(id, x_fixed[i], x[i])


def Coulomb_force(xi, xj):
    dx = xj[0] - xi[0]
    dy = xj[1] - xi[1]
    ds2 = dx * dx + dy * dy
    ds = math.sqrt(ds2)
    ds3 = ds2 * ds
    if ds3 == 0.0:
        const = 0
    else:
        const = beta / (ds2 * ds)
    return [-const * dx, -const * dy]


def Hooke_force(xi, xj, dij):
    dx = xj[0] - xi[0]
    dy = xj[1] - xi[1]
    ds = math.sqrt(dx * dx + dy * dy)
    dl = ds - dij
    const = k * dl / ds
    return [const * dx, const * dy]


def move():
    ekint = [0.0, 0.0]
    for i in xrange(m):
        Fx = 0.0
        Fy = 0.0
        for j in xrange(m):
            if i == j:
                Fij = Hooke_force(x[i], x_fixed[j], dij)
            else:
                Fij = Coulomb_force(x[i], x[j])

            Fx += Fij[0]
            Fy += Fij[1]
        v[i][0] = (v[i][0] + alpha * Fx * delta_t) * eta
        v[i][1] = (v[i][1] + alpha * Fy * delta_t) * eta
        ekint[0] = ekint[0] + alpha * (v[i][0] * v[i][0])
        ekint[1] = ekint[1] + alpha * (v[i][1] * v[i][1])

    print "total kinetic energy: %lf" % math.sqrt(ekint[0] * ekint[0] + ekint[1] * ekint[1])

    for i in xrange(m):
        x[i][0] += v[i][0] * delta_t
        x[i][1] += v[i][1] * delta_t
        move_oval(ids[i], xi=x[i])

    for i in xrange(m):
        id = lids[i]
        move_line(id, x_fixed[i], x[i])

    root.after(UPDATE_FREQ, move)

root.after(UPDATE_FREQ, move)

root.mainloop()
