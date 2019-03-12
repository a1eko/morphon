from __future__ import division, absolute_import
from builtins import sum, list, max, filter, dict, abs
from morphon import Node, rotation
from .morph import _angle_between
import numpy as np
import math


def _norm(v):
    return math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])


def _sample(c, r, h):
    link = list(_norm(c0-c1) for (c0, c1) in zip(c[:-1], c[1:]))
    path = list(sum(link[:i]) for i in range(len(link)+1))
    length = path[-1]
    n = int(length / h)
    s = list(np.linspace(0, n*h, n, endpoint=False))
    ct = c.transpose()
    xnew = np.interp(s, path, ct[0])
    ynew = np.interp(s, path, ct[1])
    znew = np.interp(s, path, ct[2])
    cnew = np.array(list((x, y, z) for (x, y, z) in zip(xnew, ynew, znew)))
    rnew = np.interp(s, path, r)
    return cnew[1:], rnew[1:]


def resample(morph, res):
    m = morph.copy()
    count = max(m.nodes.keys())+1
    for stem in m.stems():
        for sec in m.sections(stem):
            length = sum(m.length(item) for item in sec)
            if length > res:
                pttype = m.type(sec[0])
                parent = m.parent(sec[0])
                idents = [parent] + sec
                sec_coords = np.array(list(m.coord(item) for item in idents))
                sec_radii = np.array(list(m.radius(item) for item in idents))
                new_coords, new_radii = _sample(sec_coords, sec_radii, res)
                if len(new_coords) > 0:
                    for item in sec[:-1]:
                        m.remove(item)
                    item = sec[-1]
                    for coord, radius in zip(new_coords[-1::-1], new_radii[-1::-1]):
                        node = Node(count)
                        node.value = [pttype, coord, radius]
                        m.insert(item, node)
                        item = count
                        count += 1
                else:
                    if len(sec) > 1:
                        for item in sec[:-1]:
                            m.remove(item)
    m.renumber()
    return m


def select(m, secpts, types=[], orders=[], degrees=[], sec=False):
    for points in secpts:
        if types:
            if sec:
                points = points if m.type(points[-1]) in types else []
            else:
                points = list(filter(lambda i: m.type(i) in types, points))
        if orders:
            if sec and len(points) > 1:
                points = points if m.order(points[0]) in orders else []
            else:
                points = list(filter(lambda i: m.order(i) in orders, points))
        if degrees:
            if sec and len(points) > 1:
                points = points if m.degree(points[0]) in degrees else []
            else:
                points = list(filter(lambda i: m.degree(i) in degrees, points))
        if points:
            yield points


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def plot(m=None, secpts=[], ax=None, fig=None, **kwargs):
    if ax is None:
        if fig is None:
            fig = plt.figure(figsize=(8, 8))
        ax = Axes3D(fig)
        ax.xaxis.pane.set_edgecolor('w')
        ax.yaxis.pane.set_edgecolor('w')
        ax.zaxis.pane.set_edgecolor('w')
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.set_proj_type('ortho')
        ax.set_xlabel('x')
        ax.set_zlabel('y')
        ax.set_ylabel('z')
        ax.grid(False)
    else:
        fig = ax.figure
    for sec in secpts:
        c = np.array(list(m.coord(item) for item in sec))
        x, y, z = c[:,0], c[:,1], c[:,2]
        ax.plot(x, z, y, **kwargs)
    xmin, ymin, zmin = ax.xy_dataLim.xmin, ax.xy_dataLim.ymin, ax.zz_dataLim.xmin
    xmax, ymax, zmax = ax.xy_dataLim.xmax, ax.xy_dataLim.ymax, ax.zz_dataLim.xmax
    smax = max(max(ax.xy_dataLim.size), max(ax.zz_dataLim.size))
    ax.set_xlim((xmin+xmax-smax)/2, (xmin+xmax+smax)/2)
    ax.set_zlim((ymin+ymax-smax)/2, (ymin+ymax+smax)/2)
    ax.set_ylim((zmin+zmax-smax)/2, (zmin+zmax+smax)/2)
    return fig, ax


def repair_cut(m, ident, src, flip_z=False):
    root = m.root()
    n = m.copy(src)
    c0 = n.coord(src)
    cs = np.array(list(m.coord(item) for item in m.section(src)))
    c1 = np.mean(cs.transpose(), axis=1)
    v0 = c1 - c0
    if flip_z:
        n.translate(-c0)
        n.rotate(v0, math.pi)
        n.translate(c0)
    r0 = m.radius(src)
    r1 = m.radius(ident)
    n.scale(r1/r0, diam=True)
    c0 = m.coord(list(m.section(ident, reverse=True))[-1])
    c1 = m.coord(ident)
    v1 = c1 - c0
    axis, angle = rotation(v0, v1)
    n.rotate(axis, angle)
    c0 = n.coord(n.root())
    c1 = m.coord(ident)
    c1 += c1 - m.coord(m.parent(ident))
    n.translate(c1-c0)
    m.graft(ident, n)


def stretch(m, section, factor, preserve_length=True):
    if len(section) > 1:
        cs = np.array(list(m.coord(ident) for ident in section))
        c1 = np.mean(cs.transpose(), axis=1)
        c0 = m.coord(section[0])
        direction = c1 - c0
        #direction /= np.linalg.norm(direction)
        direction /= _norm(direction)
        direction *= abs(factor)
        for ident in section[1:]:
            parent = m.parent(ident)
            v = m.coord(ident) - m.coord(parent)
            #vnorm = np.linalg.norm(v)
            vnorm = _norm(v)
            vnew = v + direction*vnorm
            if preserve_length:
                vleng = np.linalg.norm(vnew)
                if vleng > 1e-6:
                    vnew /= vleng
                vnew *= vnorm
            m.translate(vnew-v, ident=ident)


def _contains(a, b):
    return any(map(lambda x: x in b, a))


def _angle(m, ident):
    theta = 0.0
    if m.is_bifurcation(ident):
        ch1, ch2 = m.children(ident)
        v1 = m.coord(ch1) - m.coord(ident)
        v2 = m.coord(ch2) - m.coord(ident)
        theta = _angle_between(v1, v2)
    return theta/math.pi*180


std_features = ['angle', 'area', 'asym', 'degree', 'depth', 'diam', 'dist',
    'height', 'length', 'nbifs', 'npoints', 'nsecs', 'ntips', 'order', 'part',
    'path', 'seclen', 'tort', 'volume', 'width']


def meter(m, ident=None, features=[]):
    root = m.root()
    ident = ident if ident is not None else root
    if not features:
        features = std_features
    c0 = m.coord(root)
    nodes = dict()
    tips = list()
    bifs = list()
    secs = list()
    part = list()
    asym = list()
    diam = list()
    angle = list()
    cmin = np.finfo(1.0).max * np.ones(3)
    cmax = np.finfo(1.0).min * np.ones(3)
    rmax = 0.0
    seclen = 0.0
    npoints = 0
    c1 = m.coord(ident)
    #
    if _contains(features, ['width', 'height', 'depth', 'dist', 'path', 'order',
        'seclen', 'tort', 'nbifs', 'nsecs', 'npoints', 'diam']):
        for item in m.traverse(ident):
            if _contains(features, ['npoints']):
                npoints += 1
            if _contains(features, ['diam']):
                diam.append(m.diam(item))
            c = m.coord(item)
            if _contains(features, ['width', 'height', 'depth']):
                cmin = np.minimum(c, cmin)
                cmax = np.maximum(c, cmax)
            if _contains(features, ['dist']):
                r = _norm(c - c0)
                rmax = max(r, rmax)
            nodes[item] = dict()
            parent = m.parent(item)
            if parent in nodes:
                seglen = m.length(item)
                if _contains(features, ['order']):
                    order = nodes[parent]['order']
                if _contains(features, ['path']):
                    path = nodes[parent]['path'] + seglen
                seclen += seglen
            else:
                if _contains(features, ['path']):
                    path = m.distance(item)
                if _contains(features, ['order']):
                    order = m.order(item)
            if _contains(features, ['path']):
                nodes[item]['path'] = path
            if m.is_bifurcation(item):
                bifs.append(item)
                if _contains(features, ['order']):
                    order = order + 1 if parent in nodes else order
            if c1 is None:
                c1 = m.coord(parent)
            if m.is_bifurcation(item) or m.is_leaf(item):
                if _contains(features, ['seclen', 'tort', 'nsecs']):
                    chord = _norm(c - c1)
                    if seclen == 0.0:
                        seclen = m.length(item)
                    secs.append(dict(length=seclen, tort=chord/seclen))
                    seclen = 0.0
                c1 = None
            if _contains(features, ['order']):
                nodes[item]['order'] = order
            if m.is_leaf(item):
                tips.append(item)
    #
    if _contains(features, ['length', 'volume', 'area', 'degree', 'part', 'asym']):
        if not tips:
            tips = list(m.leaves(ident))
            for item in m.traverse(ident):
                nodes[item] = dict()
        for tip in tips:
            for item in m.traverse(tip, reverse=True):
                if not m.children(item):
                    if _contains(features, ['length', 'asym']):
                        nodes[item]['length'] = 0.0
                    if _contains(features, ['volume']):
                        nodes[item]['volume'] = 0.0
                    if _contains(features, ['area']):
                        nodes[item]['area'] = 0.0
                    if _contains(features, ['degree', 'part']):
                        nodes[item]['degree'] = 1
                else:
                    degree = 0
                    length = 0.0
                    volume = 0.0
                    area = 0.0
                    for child in m.children(item):
                        if 'length' in nodes[child]:
                            if _contains(features, ['length', 'asym']):
                                length += nodes[child]['length']
                        if 'volume' in nodes[child]:
                            if _contains(features, ['volume']):
                                volume += nodes[child]['volume']
                        if 'area' in nodes[child]:
                            if _contains(features, ['area']):
                                area += nodes[child]['area']
                        if 'degree' in nodes[child]:
                            if _contains(features, ['degree', 'part']):
                                degree += nodes[child]['degree']
                        if _contains(features, ['length', 'asym']):
                            length += m.length(child)
                        if _contains(features, ['volume']):
                            volume += m.volume(child)
                        if _contains(features, ['area']):
                            area += m.area(child)
                    if _contains(features, ['length', 'asym']):
                        nodes[item]['length'] = length
                    if _contains(features, ['volume']):
                        nodes[item]['volume'] = volume
                    if _contains(features, ['area']):
                        nodes[item]['area'] = area
                    if _contains(features, ['degree', 'part']):
                        nodes[item]['degree'] = degree
                if item == ident:
                    break
    #
    if _contains(features, ['part', 'asym', 'angle']):
        if not bifs:
            bifs = list(m.bifurcations(ident))
        for bif in bifs:
            ch1, ch2 = m.children(bif)
            if _contains(features, ['part']):
                s1 = nodes[ch1]['degree']
                s2 = nodes[ch2]['degree']
                pa = abs(s1-s2)/(s1+s2-2) if s1+s2>2 else 0.0
                part.append(pa)
            if _contains(features, ['asym']):
                s1 = nodes[ch1]['length'] + m.length(ch1)
                s2 = nodes[ch2]['length'] + m.length(ch2)
                ds = abs(s1-s2)
                ma = (ds/s2 if s1<s2 else ds/s1) if s1 != s2 else 0.0
                asym.append(ma)
            if _contains(features, ['angle']):
                angle.append(_angle(m, bif))
    #
    mm = dict()
    if _contains(features, ['width', 'height', 'depth']):
        mm['bounds'] = cmin, cmax
    if 'width' in features:
        mm['width'] = (cmax-cmin)[0]
    if 'height' in features:
        mm['height'] = (cmax-cmin)[1]
    if 'depth' in features:
        mm['depth'] = (cmax-cmin)[2]
    if 'dist' in features:
        mm['dist'] = rmax
    if 'path' in features:
        mm['path'] = max(list(nodes[item]['path'] for item in nodes))
    if 'order' in features:
        mm['order'] = max(list(nodes[item]['order'] for item in nodes))
    if 'degree' in features:
        mm['degree'] = max(list(nodes[item]['degree'] for item in nodes))
    if _contains(features, ['seclen', 'nsecs']):
        mm['seclen'] = list(x['length'] for x in secs)
    if 'tort' in features:
        mm['tort'] = list(x['tort'] for x in secs)
    if 'nbifs' in features:
        mm['nbifs'] = len(bifs)
    if 'nsecs' in features:
        mm['nsecs'] = len(secs)
    if 'ntips' in features:
        mm['ntips'] = len(tips)
    if 'length' in features:
        mm['length'] = nodes[ident]['length']
    if 'volume' in features:
        mm['volume'] = nodes[ident]['volume']
    if 'area' in features:
        mm['area'] = nodes[ident]['area']
    if 'part' in features:
        mm['part'] = part
    if 'asym' in features:
        mm['asym'] = asym
    if 'angle' in features:
        mm['angle'] = angle
    if 'npoints' in features:
        mm['npoints'] = npoints
    if 'diam' in features:
        mm['diam'] = diam
    return mm


def _sholl_crossings(m, i, h):
    p = m.parent(i)
    r1 = m.distance(p, radial=True)
    r2 = m.distance(i, radial=True)
    k1 = int(r1/h)
    k2 = int(r2/h)
    return k2-k1, k1


def sholl(m, idents, step=10):
    rmax = max(m.distance(i, radial=True) for i in idents)
    radx = np.array([k*step for k in range(int(rmax/step))])
    crox = np.zeros(int(rmax/step), dtype=int)
    for ident in idents:
        ncross, icross = _sholl_crossings(m, ident, step)
        if ncross > 0:
            for k in range(ncross):
                crox[icross+k] += 1
        elif ncross < 0:
            for k in range(1,-ncross):
                crox[icross-k] += 1
    return radx, crox
