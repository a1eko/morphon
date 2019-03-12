from __future__ import division, absolute_import
from builtins import list, dict, filter, all, sum, min, max, range
import numpy as np
import math
import copy

from .nodes import Node, Nodes
#from .lib import _norm
from . import swc


def _norm(v):
    return math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])


def _rotation_matrix(axis, angle):
    axis = np.asarray(axis)
    angle = np.asarray(angle)
    #axis = axis/np.sqrt(np.dot(axis, axis)) if np.dot(axis, axis) > 1e-6 else axis
    axis = axis/math.sqrt(np.dot(axis, axis)) if np.dot(axis, axis) > 1e-6 else axis
    a = math.cos(angle/2)
    #b, c, d = -axis*np.sin(angle/2)
    b, c, d = -axis*math.sin(angle/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def _unit_vector(vector):
    norm = _norm(vector)
    return vector / norm if norm > 1e-6 else vector # FIXME


def _angle_between(v1, v2):
    v1_u = _unit_vector(v1)
    v2_u = _unit_vector(v2)
    #angle = np.arccos(np.dot(v1_u, v2_u))
    angle = math.acos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return math.pi
    return angle


def rotation(v1, v2):
    angle = _angle_between(v1, v2)
    axis = _unit_vector(np.cross(v1, v2)) if angle != 0.0 and angle != math.pi else np.array(v1)*0
    return axis, angle


class Point(object):
    UNDEF, SOMA, AXON, DEND, APIC = range(5)
    TYPES = (SOMA, AXON, DEND, APIC)
    NEURITES = (AXON, DEND, APIC)


class Morph(Nodes):
    def __init__(self):
        Nodes.__init__(self)

    def load(self, source):
        swc.load(self, source)

    def save(self, target, header=''):
        swc.save(self, target, header=header)

    def type(self, ident):
        t, xyz, r = self.value(ident)
        return t

    def coord(self, ident):
        t, xyz, r = self.value(ident)
        return copy.deepcopy(xyz)

    def radius(self, ident):
        t, xyz, r = self.value(ident)
        return r

    def diam(self, ident):
        t, xyz, r = self.value(ident)
        return 2*r

    def length(self, ident):
        parent = self.parent(ident)
        if parent is not None:
            c1 = self.coord(ident)
            c0 = self.coord(parent)
            L = _norm(c0-c1)
        else:
            L = 0
        return L

    def area(self, ident):
        parent = self.parent(ident)
        if parent is not None:
            h = self.length(ident)
            r1 = self.radius(ident)
            r0 = self.radius(parent)
            A = math.pi*(r0+r1)*math.sqrt((r0-r1)*(r0-r1) + h*h)
        else:
            A = 0
        return A

    def volume(self, ident):
        parent = self.parent(ident)
        if parent is not None:
            h = self.length(ident)
            r1 = self.radius(ident)
            r0 = self.radius(parent)
            V = math.pi/3*(r0*r0 + r0*r1 + r1*r1)*h
        else:
            V = 0
        return V

    def is_soma(self, ident):
        return self.type(ident) == Point.SOMA

    def is_axon(self, ident):
        return self.type(ident) == Point.AXON

    def is_dend(self, ident):
        return self.type(ident) == Point.DEND

    def is_apic(self, ident):
        return self.type(ident) == Point.APIC

    def stems(self):
        return filter(lambda i: not self.is_soma(i), self.children(self.root()))

    def distance(self, ident, radial=False, reverse=True):
        if radial:
            c0 = self.coord(self.root())
            c1 = self.coord(ident)
            dist = _norm(c0-c1)
        else:
            dist = sum(self.length(item) for item in self.traverse(ident, reverse=reverse))
        return dist

    def move(self, ident, shift):
        self.nodes[ident].value[1] += np.asarray(shift)

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

    def copy(self, ident=None):
        if ident:
            m = Morph()
            for item in self.traverse(ident):
                m.nodes[item] = copy.deepcopy(self.nodes[item])
            m.nodes[ident].parent = None
        else:
            m = copy.deepcopy(self)
        return m

    def check(self):
        rep = Nodes.check(self)
        rep['increasing_id_order'] = all(self.parent(item) < item
            for item in filter(lambda i: not self.is_root(i), self.nodes))
        rep['valid_types'] = all(self.type(item) in Point.TYPES for item in self.nodes)
        rep['nonzero_diam'] = all(self.diam(item) != 0 for item in self.nodes)
        rep['sequential_ids'] = sorted(self.nodes) == list(range(min(self.nodes), max(self.nodes)+1))
        soma = list(filter(self.is_soma, self.nodes))
        if soma:
            rep['unit_id_in_soma'] = 1 in soma
            rep['unit_id_is_root'] = 1 == self.root()
            path = list()
            for ident in self.stems():
                path.extend(self.traverse(ident))
            rep['unit_id_is_origin'] = set(path) == set(self.nodes).difference(soma)
        return rep
