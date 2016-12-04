from __future__ import print_function, division
from morphon import Morpho
from math import pi

import matplotlib.pyplot as plt
import numpy as np


title = 'Neuron plots'
morphology  = 'in.swc'
plot = 'out3b.pdf'

def colorize(neurite):
    if neurite == 'soma':
        c = 'black'
    elif neurite == 'axon':
        c = 'red'
    elif neurite == 'dend':
        c = 'green'
    elif neurite == 'apic':
        c = 'blue'
    else:
        c = 'cyan'  # unknown
    return c

def plot3d(m, ax):
    ax.set_axis_off()
    for branch in m.branches():
        cx = [m.coord(i)[0] for i in branch]
        cy = [m.coord(i)[1] for i in branch]
        cz = [m.coord(i)[2] for i in branch]
        parent = m.parent(branch[0])
        if parent is not None:
            cx.insert(0, m.coord(parent)[0])
            cy.insert(0, m.coord(parent)[1])
            cz.insert(0, m.coord(parent)[2])
        cx = np.array(cx)
        cy = np.array(cy)
        cz = np.array(cz)
        x = cx+0.5*cz
        y = cy+0.5*cz
        neurite = m.neurite(branch[0])
        c = colorize(neurite)
        lw = 5 if neurite is 'soma' else 1
        ax.plot(x, y, color=c, linewidth=lw)

def plot2d(m, ax, projection, xlim=None, ylim=None):
    ax.set_axis_on()
    for branch in m.branches():
        cx = [m.coord(i)[0] for i in branch]
        cy = [m.coord(i)[1] for i in branch]
        cz = [m.coord(i)[2] for i in branch]
        parent = m.parent(branch[0])
        if parent is not None:
            cx.insert(0, m.coord(parent)[0])
            cy.insert(0, m.coord(parent)[1])
            cz.insert(0, m.coord(parent)[2])
        cx = np.array(cx)
        cy = np.array(cy)
        cz = np.array(cz)
        if projection is 'xy':
            x = cx
            y = cy
            ax.set_xlabel('x')
            ax.set_ylabel('y')
        elif projection is 'yz':
            x = cz
            y = cy
            ax.set_xlabel('z')
            ax.set_ylabel('y')
        else:
            x = cx
            y = cz
            ax.set_xlabel('x')
            ax.set_ylabel('z')
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        neurite = m.neurite(branch[0])
        c = colorize(neurite)
        lw = 5 if neurite is 'soma' else 1
        ax.plot(x, y, color=c, linewidth=lw)

print(title)
print('neuron:', morphology)

m = Morpho(morphology)
(xmin, ymin, zmin), (xmax, ymax, zmax) = m.bounds()
sx, sy, sz = xmax-xmin, ymax-ymin, zmax-zmin
smax = max([sx, sy, sz])

fig = plt.figure(figsize=(8,8))
fig.suptitle('Neuron ' + morphology)

ax = plt.subplot(2, 2, 1)
plt.plot([xmin+zmax/2, xmin+zmax/2], [ymin+zmax/2, ymax+zmax/2], c='k')
plt.plot([xmin+zmax/2, xmax+zmax/2], [ymin+zmax/2, ymin+zmax/2], c='k')
plt.plot([xmin+zmax/2, xmin], [ymin+zmax/2, ymin], c='k')
plot3d(m, ax)

ax = plt.subplot(2, 2, 2)
ax.set_xlim(((xmin+xmax-smax)/2, (xmin+xmax+smax)/2))
ax.set_ylim(((ymin+ymax-smax)/2, (ymin+ymax+smax)/2))
plot2d(m, ax, projection='xy')

ax = plt.subplot(2, 2, 3)
ax.set_xlim(((zmin+zmax-smax)/2, (zmin+zmax+smax)/2))
ax.set_ylim(((ymin+ymax-smax)/2, (ymin+ymax+smax)/2))
plot2d(m, ax, projection='yz')

ax = plt.subplot(2, 2, 4)
ax.set_xlim(((xmin+xmax-smax)/2, (xmin+xmax+smax)/2))
ax.set_ylim(((zmin+zmax-smax)/2, (zmin+zmax+smax)/2))
plot2d(m, ax, projection='xz')

plt.tight_layout()
plt.savefig(plot)
print('image saved to', plot)
