from __future__ import division
from morphon import Morpho
from math import pi
import numpy as np

def measure(m):
    morphometrics = {}
    all_nodes = [i for i in m.traverse()]
    all_branches = [b for b in m.branches()]
    neurites = set(m.neurite(i) for i in all_nodes)
    for neurite in neurites:
        measures = {}
        branches = filter(lambda b: m.neurite(b[0])==neurite, all_branches)
        nodes = filter(lambda i: m.neurite(i)==neurite, all_nodes)
        bifurcations = filter(lambda i: m.is_fork(i), nodes)
        tips = filter(lambda i: m.is_leaf(i), nodes)
        stems = filter(lambda b: m.order(b[0])==1, branches)
        diams =[m.diam(i) for i in nodes]
        area = sum(m.area(i) for i in nodes)
        length = sum(m.length(i) for i in nodes)
        volume = sum(m.volume(i) for i in nodes)
        measures['area'] = area
        measures['length'] = length
        measures['volume'] = volume
        measures['extent euclidean'] = [s for s in m.size(idents=nodes)]
        measures['extent radial'] = max(m.distance(i, radial=True) for i in nodes)
        measures['extent path'] = max(m.distance(i) for i in tips)
        measures['diameters'] = np.mean(diams), np.std(diams), np.min(diams), np.max(diams)
        measures['diameter effective'] = area/(pi*length)
        measures['order'] = max(m.order(i) for i in tips)
        measures['stems'] = len(stems)
        measures['bifurcations'] = len(bifurcations)
        measures['tips'] = len(tips)
        measures['branches'] = len(branches)
        morphometrics[neurite] = measures
    return morphometrics
    

def _crossings(m, i, h):
    p = m.parent(i)
    r1 = m.distance(p, radial=True)
    r2 = m.distance(i, radial=True)
    k0 = int(r1/h)
    kn = int(r2/h)
    return kn-k0, k0


def scholl(m, neurite, h=10):
    nodes = filter(lambda i: m.neurite(i)==neurite, m.traverse())
    rmax = max([m.distance(i, radial=True) for i in nodes])
    radx = np.array([k*h for k in range(int(rmax/h))])
    crox = np.zeros(int(rmax/h), dtype=int)
    for i in nodes:
        ncross, icross = _crossings(m, i, h)
        if ncross > 0:
            for k in range(ncross):
                crox[icross+k] += 1
        elif ncross < 0:
            for k in range(-ncross):
                crox[icross-k] += 1
    return radx, crox


def plot(m, ax, projection='xy', neurites=[], orders=[], degrees=[], color='b', linewidth=1, equal_scales=False):
    for section in m.sections(with_parent=True, neurites=neurites, orders=orders, degrees=degrees):
        x, y, z = m.coords(section)
	if linewidth == 'diam':
	    lw = (m.diam(section[-len(section)+1]) + m.diam(section[-1])) * 0.5
	else:
	    lw = linewidth
	if equal_scales:
	    (xmin, ymin, zmin), (xmax, ymax, zmax) = m.bounds()
	    sx, sy, sz = xmax-xmin, ymax-ymin, zmax-zmin
	    smax = max([sx, sy, sz])
	if projection == 'xy':
	    if equal_scales:
	        ax.set_xlim(((xmin+xmax-smax)/2, (xmin+xmax+smax)/2))
		ax.set_ylim(((ymin+ymax-smax)/2, (ymin+ymax+smax)/2))
	    ax.plot(x, y, color=color, linewidth=lw)
	elif projection == 'yz':
	    if equal_scales:
	        ax.set_xlim(((ymin+ymax-smax)/2, (ymin+ymax+smax)/2))
		ax.set_ylim(((zmin+zmax-smax)/2, (zmin+zmax+smax)/2))
	    ax.plot(y, z, color=color, linewidth=lw)
	elif projection == 'xz':
	    if equal_scales:
	        ax.set_xlim(((xmin+xmax-smax)/2, (xmin+xmax+smax)/2))
		ax.set_ylim(((zmin+zmax-smax)/2, (zmin+zmax+smax)/2))
	    ax.plot(x, z, color=color, linewidth=lw)
	elif projection == '3d' or projection == 'xzy':
	    if equal_scales:
	        ax.set_xlim(((xmin+xmax-smax)/2, (xmin+xmax+smax)/2))
	        ax.set_zlim(((ymin+ymax-smax)/2, (ymin+ymax+smax)/2))
		ax.set_ylim(((zmin+zmax-smax)/2, (zmin+zmax+smax)/2))
	    ax.plot(x, z, y, color=color, linewidth=lw)
	elif projection == 'xyz':
	    ax.plot(x, y, z, color=color, linewidth=lw)
        else:
	    raise Exception('unknown projection ' + projection)
