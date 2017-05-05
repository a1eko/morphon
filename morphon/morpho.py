from __future__ import division

import math
import copy
import numpy as np

from tree import Tree


def _rotation_matrix(axis, angle):
    axis = np.asarray(axis)
    angle = np.asarray(angle)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(angle/2)
    b, c, d = -axis*math.sin(angle/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def _unit_vector(vector):
    norm = np.linalg.norm(vector)
    return vector / norm

def _angle_between(v1, v2):
    v1_u = _unit_vector(v1)
    v2_u = _unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle

def rotation(v1, v2):
    angle = _angle_between(v1, v2)
    axis = _unit_vector(np.cross(v1, v2)) if angle != 0.0 and angle != np.pi else np.array(v1)*0
    return axis, angle

def elevation(v):
    v2 = np.array([v[0], v[1], 0])
    angle = _angle_between(v, v2)
    return angle if v[2] > 0 else -angle


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

    def copy(self, ident=None):
        morpho = Morpho()
        for item in self.traverse(ident):
            morpho.nodes[item] = copy.deepcopy(self.nodes[item])
        if ident is not None:
            morpho.nodes[ident].parent = None
        return morpho

    def save(self, target, header=''):
        data = []
        if target.lower().endswith('.swc'):
            for ident in self.traverse():
                (n, (x, y, z), d), p = self.value(ident), self.parent(ident)
                p = p if p is not None else -1
                r = d/2
                data.append([ident, n, x, y, z, r, p])
            np.savetxt(target, data, fmt='%d %d %g %g %g %g %d', header=header)
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
        if parent is not None:
            c1 = self.coord(ident)
            c0 = self.coord(parent)
            L = np.linalg.norm(c0-c1)
        else:
            L = self.diam(ident)
        return L

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

    def distance(self, ident, radial=False, reverse=True):
        if radial:
            c0 = self.coord(self.root())
            c1 = self.coord(ident)
            dist = np.linalg.norm(c0-c1)
        else:
            dist = sum(self.length(item) for item in self.traverse(ident, reverse=reverse))
        return dist

    def angle(self, ident):
        theta = 0.0
        node = self.nodes[ident]
        if not self.is_root(ident):
            if self.is_bifurcation(ident):
                c0 = self.coord(node.children[0])
                c1 = self.coord(node.children[1])
                v0 = c0 - self.coord(ident)
                v1 = c1 - self.coord(ident)
                theta = _angle_between(v0, v1)
            elif not self.is_leaf(ident):
                c0 = self.coord(node.children[0])
                c1 = self.coord(node.parent)
                v0 = c0 - self.coord(ident)
                v1 = self.coord(ident) - c1
                theta = _angle_between(v0, v1)
        return theta

    def curvature(self, ident):
        scalar_curvature = 0.0
        if not self.is_root(ident) and not self.is_bifurcation(ident) and not self.is_leaf(ident):
            node = self.nodes[ident]
            A = self.coord(node.parent)
            B = self.coord(ident)
            C = self.coord(node.children[0])
            a = np.linalg.norm(C-B)
            b = np.linalg.norm(C-A)
            c = np.linalg.norm(B-A)
            s = (a+b+c)/2
            q = np.sqrt(abs(s*(s-a)*(s-b)*(s-c)))
            scalar_curvature = 1 / (a*b*c/4/q) if abs(q) > 1e-9 else 0.0
        return scalar_curvature

    def increment(self, ident, axis=2):
        du = 0.0
        if axis not in [0, 1, 2]:
            if type(axis) is str:
                axis = axis.lower()
                axis = {'x': 0, 'y': 1, 'z': 2}[axis]
            else:
                raise Error('incorrect axis ' + axis)
        parent = self.parent(ident)
        if parent is not None and not self.is_root(parent):
            c1 = self.coord(ident)
            c0 = self.coord(parent)
            du = c1[axis] - c0[axis]
        return du

    def rel_increment(self, ident, axis=2):
        dl = 1.0
        du = self.increment(ident, axis=axis)
        parent = self.parent(ident)
        if parent is not None:
            dl = self.length(parent)
        return du / dl

    def is_jump(self, ident, axis=2, increment_thresh=5, rel_increment_thresh=3):
        jump = False
        if not self.is_leaf(ident):
            child = self.nodes[ident].children[0]
            if (abs(self.increment(ident, axis=axis)) > increment_thresh
                and abs(self.rel_increment(ident, axis=axis)) > rel_increment_thresh
                and abs(self.rel_increment(child, axis=axis)) < rel_increment_thresh):
                jump = True
        return jump

    def jumps(self, ident=None, idents=[], axis=2, increment_thresh=5, rel_increment_thresh=3):
        if ident is None:
            ident = self.root()
        if not idents:
            idents = self.traverse(ident)
        items = [item for item in idents
            if self.is_jump(item, axis=axis, increment_thresh=increment_thresh, rel_increment_thresh=rel_increment_thresh)]
        return items

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

    def translate(self, shift, ident=None):
        for item in self.traverse(ident):
            self.nodes[item].value[1] += np.asarray(shift)

    def scale(self, factor, coord=True, diam=False, ident=None):
        for item in self.traverse(ident):
            if coord:
                self.nodes[item].value[1] *= np.asarray(factor)
            if diam:
                self.nodes[item].value[2] *= np.asarray(factor).mean()

    def rotate(self, axis, angle, ident=None):
        for item in self.traverse(ident):
            c = self.coord(item)
            self.nodes[item].value[1] = np.dot(_rotation_matrix(axis, angle), c)

    def sections(self, ident=None, reverse=False, with_parent=False,
        neurites=[], orders=[], degrees=[]):
        selected = [b for b in self.branches(ident=ident, reverse=reverse)]
        head = 0 if not with_parent else 1
        for section in selected:
            if with_parent:
                parent = self.parent(section[0])
                if parent is not None:
                    section.insert(0, parent)
                else:
                    section.insert(0, section[0])
        if neurites:
            selected= filter(lambda b: self.neurite(b[head]) in neurites, selected)
        if orders:
            selected= filter(lambda b: self.order(b[head]) in orders, selected)
        if degrees:
            selected= filter(lambda b: self.degree(b[head]) in degrees, selected)
        return selected

    def points(self, ident=None, reverse=False, neurites=[], orders=[], degrees=[]):
        selected = [i for i in self.traverse(ident=ident, reverse=reverse)]
        if neurites:
            selected= filter(lambda i: self.neurite(i) in neurites, selected)
        if orders:
            selected= filter(lambda i: self.order(i) in orders, selected)
        if degrees:
            selected= filter(lambda i: self.degree(i) in degrees, selected)
        return selected

    def tips(self, ident=None, reverse=False, neurites=[], orders=[], degrees=[]):
        idents = self.points(ident=ident, reverse=False,
            neurites=neurites, orders=orders, degrees=degrees)
        selected = filter(lambda i: self.is_leaf(i), idents)
        return selected

    def stems(self, ident=None, reverse=False, neurites=[], degrees=[]):
        stem_sections = self.sections(ident=ident, reverse=reverse, neurites=neurites, orders=[1], degrees=degrees)
        stem_sections = filter(lambda s: self.neurite(s[0]) != 'soma', stem_sections)
        return [s[0] for s in stem_sections]

    def bifurcations(self, ident=None, reverse=False, neurites=[], orders=[], degrees=[]):
        idents = self.points(ident=ident, reverse=reverse,
            neurites=neurites, orders=orders, degrees=degrees)
        selected = filter(lambda i: self.is_bifurcation(i), idents)
        return selected
