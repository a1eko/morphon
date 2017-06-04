from __future__ import division
from matplotlib.mlab import PCA
from morphon import Morpho
import numpy as np
import math


def measure(m, features=[], idents=[], ident=None, reverse=False, increment_thresh=5, rel_increment_thresh=3):
    metrics = {}
    tips = []
    stems = []
    bifurcations = []
    area = None
    length = None
    if features is None:
        features = []
    if 'number_of_stems' in features or not features:
        if not idents:
            stems = m.stems(ident, reverse=reverse)
        metrics['number_of_stems'] = len(stems) if not idents else np.nan
    if 'number_of_branches' in features or not features:
        if not idents:
            branches = [b for b in m.branches(ident, reverse=reverse) if m.neurite(b[0]) != 'soma']
        metrics['number_of_branches'] = len(branches) if not idents else np.nan
    if not idents:
        idents = [i for i in m.traverse(ident, reverse=reverse)]
    if 'pca_fracs' in features or not features:
        data = np.array([m.coord(i) for i in idents])
	results = None
	try:
            results = PCA(data)
	except:
	    pass
	metrics['pca_fracs'] = results.fracs.tolist() if results else {'pca_fracs': [np.nan, np.nan, np.nan]}
    if 'area' in features or not features:
        area = sum(m.area(i) for i in idents)
        metrics['area'] = area
    if 'length' in features or not features:
        length = sum(m.length(i) for i in idents)
        metrics['length'] = length
    if 'volume' in features or not features:
        metrics['volume'] = sum(m.volume(i) for i in idents)
    if 'root_position' in features or not features:
        root = m.root()
        metrics['root_position'] = m.coord(root).tolist() if root in idents else np.nan
    if 'center_position' in features or not features:
        positions = np.array([m.coord(i) for i in idents])
        center = sum(positions) / len(positions)
        metrics['center_position'] = center.tolist()
    if 'euclidean_extent' in features or not features:
        metrics['euclidean_extent'] = [s for s in m.size(idents=idents)]
    if 'radial_extent' in features or not features:
        metrics['radial_extent'] = max(m.distance(i, radial=True) for i in idents)
    if 'path_extent' in features or not features:
        tips = filter(lambda i: m.is_leaf(i), idents)
        metrics['path_extent'] = max(m.distance(i) for i in tips) if tips else np.nan
    if 'local_diameter' in features or not features:
        diams =[m.diam(i) for i in idents]
        metrics['local_diameter'] = np.mean(diams), np.std(diams), np.median(diams), np.min(diams), np.max(diams)
    if 'effective_diameter' in features or not features:
        if not area:
            area = sum(m.area(i) for i in idents)
        if not length:
            length = sum(m.length(i) for i in idents)
        metrics['effective_diameter'] = area/(math.pi*length)
    if 'order' in features or not features:
        orders = [m.order(i) for i in idents]
        min_order = min(orders)
        max_order = max(orders)
        metrics['order'] = min_order, max_order
    if 'degree' in features or not features:
        degrees = [m.degree(i) for i in idents]
        min_degree = min(degrees)
        max_degree = max(degrees)
        metrics['degree'] = min_degree, max_degree
    if 'bifurcation_angle' in features or not features:
        bifurcations = filter(lambda i: m.is_bifurcation(i), idents)
        if bifurcations:
            angles = [m.angle(i) for i in bifurcations]
            metrics['bifurcation_angle'] = np.mean(angles), np.std(angles), np.median(angles), np.min(angles), np.max(angles)
    if 'curvature' in features or not features:
        curvatures = [m.curvature(i) for i in idents]
        metrics['curvature'] = np.mean(curvatures), np.std(curvatures), np.median(curvatures), np.min(curvatures), np.max(curvatures)
    if 'number_of_tips' in features or not features:
        if not tips:
            tips = filter(lambda i: m.is_leaf(i), idents)
        metrics['number_of_tips'] = len(tips)
    if 'number_of_bifurcations' in features or not features:
        if not bifurcations:
            bifurcations = filter(lambda i: m.is_bifurcation(i), idents)
        metrics['number_of_bifurcations'] = len(bifurcations)
    if 'z_jumps' in features or not features:
        jumps = m.jumps(idents=idents, increment_thresh=increment_thresh, rel_increment_thresh=rel_increment_thresh)
        metrics['z_jumps'] = len(jumps)
    return metrics


def _sholl_crossings(m, i, h):
    p = m.parent(i)
    r1 = m.distance(p, radial=True)
    r2 = m.distance(i, radial=True)
    k0 = int(r1/h)
    kn = int(r2/h)
    return kn-k0, k0


def sholl(m, h=10, neurites=[], orders=[], degrees=[]):
    idents = m.points(neurites=neurites, orders=orders, degrees=degrees)
    rmax = max(m.distance(i, radial=True) for i in idents)
    radx = np.array([k*h for k in range(int(rmax/h))])
    crox = np.zeros(int(rmax/h), dtype=int)
    for i in idents:
        ncross, icross = _sholl_crossings(m, i, h)
        if ncross > 0:
            for k in range(ncross):
                crox[icross+k] += 1
        elif ncross < 0:
            for k in range(-ncross):
                crox[icross-k] += 1
    return radx, crox


def _plot_projection(ax, bounds, x, y, z=None, projection='xy', equal_scales=False, **kwargs):
    if equal_scales:
        (xmin, ymin, zmin), (xmax, ymax, zmax) = bounds
        sx, sy, sz = xmax-xmin, ymax-ymin, zmax-zmin
        smax = max([sx, sy, sz])
    if projection == 'xy':
        if equal_scales:
            ax.set_xlim(((xmin+xmax-smax)/2, (xmin+xmax+smax)/2))
            ax.set_ylim(((ymin+ymax-smax)/2, (ymin+ymax+smax)/2))
        ax.plot(x, y, **kwargs)
    elif projection == 'yz':
        if equal_scales:
            ax.set_xlim(((ymin+ymax-smax)/2, (ymin+ymax+smax)/2))
            ax.set_ylim(((zmin+zmax-smax)/2, (zmin+zmax+smax)/2))
        ax.plot(y, z, **kwargs)
    elif projection == 'xz':
        if equal_scales:
            ax.set_xlim(((xmin+xmax-smax)/2, (xmin+xmax+smax)/2))
            ax.set_ylim(((zmin+zmax-smax)/2, (zmin+zmax+smax)/2))
        ax.plot(x, z, **kwargs)
    elif projection == '3d' or projection == 'xzy':
        if equal_scales:
            ax.set_xlim(((xmin+xmax-smax)/2, (xmin+xmax+smax)/2))
            ax.set_zlim(((ymin+ymax-smax)/2, (ymin+ymax+smax)/2))
            ax.set_ylim(((zmin+zmax-smax)/2, (zmin+zmax+smax)/2))
        ax.plot(x, z, y, **kwargs)
    elif projection == 'xyz':
        ax.plot(x, y, z, **kwargs)
    else:
        raise Exception('unknown projection ' + projection)


def plot(m, ax, projection='xy', neurites=[], orders=[], degrees=[], idents=[], equal_scales=False, **kwargs):
    if not idents:
        for section in m.sections(with_parent=True, neurites=neurites, orders=orders, degrees=degrees):
            x, y, z = m.coords(section)
            _plot_projection(ax, m.bounds(), x, y, z=z, projection=projection, equal_scales=equal_scales, **kwargs)
    else:
        x, y, z = m.coords(idents)
        _plot_projection(ax, m.bounds(), x, y, z=z, projection=projection, equal_scales=equal_scales, **kwargs)


def correct_jumps(m, jumps, axis=2):
    if axis not in [0, 1, 2]:
        if type(axis) is str:
            axis = axis.lower()
            axis = {'x': 0, 'y': 1, 'z': 2}[axis]
        else:
            raise Error('incorrect axis ' + axis)
    for ident in jumps:
        parent = m.parent(ident)
        c1 = m.coord(ident)
        c0 = m.coord(parent)
        shift = c0[axis]-c1[axis]
        corr = [0, 0, 0]
        corr[axis] = shift
        m.translate(corr, ident=ident)


def tortuosity(m, section):
    c0 = m.coord(section[0])
    c1 = m.coord(section[-1])
    dist = np.linalg.norm(c1-c0)
    leng = sum(m.length(ident) for ident in section)
    return dist / leng if abs(leng) > 1e-6 else 1.0


def stretch(m, section, scaling, preserve_length=True, preserve_bifurcation_angle=True):
    start = 2 if preserve_bifurcation_angle else 1
    if len(section) > start:
        positions = np.array([m.coord(ident) for ident in section])
        center = sum(positions) / len(positions)
        origin = m.coord(section[0])
        direction = center - origin
        direction /= np.linalg.norm(direction)
        direction *= abs(scaling)
        for ident in section[start:]:
            parent = m.parent(ident)
            v = m.coord(ident) - m.coord(parent)
            vnorm = np.linalg.norm(v)
            vnew = v + direction*vnorm
            if preserve_length:
                vnew /= np.linalg.norm(vnew)
                vnew *= vnorm
            shift = vnew - v
            m.translate(shift, ident=ident)
        ident = section[-1]
        if m.is_bifurcation(ident):
            c0 = m.coord(ident)
            for child in m.nodes[ident].children:
                c1 = m.coord(child)
                m.translate(c0-c1, ident=child)
