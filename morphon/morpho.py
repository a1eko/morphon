from __future__ import division

import math
import copy
import numpy as np

from tree import Tree


def rotation_matrix(axis, theta):
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2)
    b, c, d = -axis*math.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

Neurite = {'soma': 1, 'axon': 2, 'dend': 3, 'apic': 4}
Neurites = {v: k for k, v in Neurite.items()}


class Error(StandardError):
    pass


class Morpho(Tree):
    def __init__(self, source=None):
        Tree.__init__(self)
        if source:
            self.load(source)

    def load(self, source):
        if source.lower().endswith('.swc'):
            data = np.loadtxt(source, dtype={'names': ('i', 'n', 'x', 'y', 'z', 'r', 'p'),
                'formats': (int, int, float, float, float, float, int)})
            for i, n, x, y, z, r, p in data:
                n = n if n in Neurites else Neurite[self.neurite(p)]
                Tree.add(self, ident=i, value=[n, np.array([x, y, z]), 2*r], parent=p)
        else:
            raise Error('unknown file format of ' + source)

    def save(self, target, header=''):
        data = []
        if target.lower().endswith('.swc'):
            for ident in self.traverse():
                (n, (x, y, z), d), p = self.value(ident), self.parent(ident)
                p = p if p is not None else -1
                r = d/2
                data.append([ident, n, x, y, z, r, p])
            np.savetxt(target, data, fmt='%d %d %.3g %.3g %.3g %.3g %d', header=header)
        else:
            raise Error('unknown file format of ' + target)

    def add(self, neurite, coord, diam, parent=None):
        ident = max(self.nodes.keys())+1 if self.nodes else None
        return Tree.add(self, ident=ident, value=[Neurite[neurite], np.array(coord), diam], parent=parent)


    def neurite(self, ident):
        n, xyz, d = self.value(ident)
        return Neurites[n]

    def coord(self, ident):
        n, xyz, d = self.value(ident)
        return copy.deepcopy(xyz)

    def coords(self, idents):
        x = [self.coord(i)[0] for i in idents]
	y = [self.coord(i)[1] for i in idents]
	z = [self.coord(i)[2] for i in idents]
        return x, y, z

    def diam(self, ident):
        n, xyz, d = self.value(ident)
        return d

    def length(self, ident):
        parent = self.parent(ident)
        c1 = self.coord(ident)
        c0 = self.coord(parent) if parent is not None else c1
        return np.linalg.norm(c0-c1)

    def area(self, ident):
        parent = self.parent(ident)
        h = self.length(ident)
        r1 = self.diam(ident)/2
        r0 = self.diam(parent)/2 if parent is not None else r1
        return math.pi*(r0+r1)*math.sqrt((r0-r1)*(r0-r1) + h*h)

    def volume(self, ident):
        parent = self.parent(ident)
        h = self.length(ident)
        r1 = self.diam(ident)/2
        r0 = self.diam(parent)/2 if parent is not None else 0
        return math.pi/3*(r0*r0 + r0*r1 + r1*r1)*h

    def distance(self, ident, radial=False):
	if radial:
            c0 = self.coord(self.root())
            c1 = self.coord(ident)
            dist = np.linalg.norm(c0-c1)
	else:
            dist = sum(self.length(item) for item in self.traverse(ident, reverse=True))
        return dist

    def bounds(self, ident=None, reverse=False, idents=[]):
	if not idents:
            if ident is None:
                ident = self.root()
            cmin = self.coord(ident)
            cmax = cmin
            for item in self.traverse(ident, reverse=reverse):
                c = self.coord(item)
                cmin = np.minimum(c, cmin)
                cmax = np.maximum(c, cmax)
	else:
            cmin = self.coord(idents[0])
            cmax = cmin
            for item in idents:
                c = self.coord(item)
                cmin = np.minimum(c, cmin)
                cmax = np.maximum(c, cmax)
        return cmin, cmax

    def size(self, ident=None, reverse=False, idents=[]):
        cmin, cmax = self.bounds(ident, reverse=reverse, idents=idents)
        return cmax-cmin

    def translate(self, ident, shift):
        self.nodes[ident].value[1] += np.asarray(shift)

    def scale(self, ident, factor, coord=True, diam=False):
        if coord:
            self.nodes[ident].value[1] *= np.asarray(factor)
        if diam:
            self.nodes[ident].value[2] *= np.asarray(factor).mean()

    def rotate(self, ident, axis, angle):
        c = self.coord(ident)
        self.nodes[ident].value[1] = np.dot(rotation_matrix(axis, angle), c)

    def sections(self, ident=None, depth=True, reverse=False, with_parent=False, 
        neurites=[], orders=[], degrees=[]):
        selected = [b for b in self.branches(ident=ident, depth=depth, reverse=reverse)]
	for section in selected:
	    if with_parent:
	        parent = self.parent(section[0]) 
		if parent is not None:
	            section.insert(0, parent)
	if neurites:
	    selected= filter(lambda b: self.neurite(b[-1]) in neurites, selected)
	if orders:
	    selected= filter(lambda b: self.order(b[-1]) in orders, selected)
	if degrees:
	    selected= filter(lambda b: self.degree(b[-1]) in degrees, selected)
        return selected
