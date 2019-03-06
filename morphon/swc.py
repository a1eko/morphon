from __future__ import absolute_import
from builtins import list
import numpy as np

from .nodes import Node, Nodes


class Error(Exception):
    pass


def _parent(p):
    if p == -1:
        p = None
    elif p is None:
        p = -1
    else:
        p = int(p)
    return p


def load(m, source):
    if source.lower().endswith('.swc'):
        nam = 'i', 't', 'x', 'y', 'z', 'r', 'p'
        fmt = 'i', 'i', 'f', 'f', 'f', 'f', 'i'
        data = np.loadtxt(source, dtype={'names': nam, 'formats': fmt})
        for i, t, x, y, z, r, p in data:
            m.add(Node(int(i), value=[int(t), np.array([x, y, z]), r], parent=_parent(p)))
        for item in m.nodes:
            m.connect(item, m.parent(item))
    else:
        raise Error('unknown file format of ' + source)


def save(m, target, header=''):
    if target.lower().endswith('.swc'):
        data = list()
        for item in m.traverse():
            (t, (x, y, z), r), p = m.value(item), m.parent(item)
            data.append([item, t, x, y, z, r, _parent(p)])
        np.savetxt(target, data, fmt='%d %d %g %g %g %g %d', header=header)
    else:
        raise Error('unknown file format of ' + target)
