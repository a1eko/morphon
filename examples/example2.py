from __future__ import print_function, division
from morphon import Morpho
from math import pi
import numpy as np
import json


title = 'Basic morphometrics'
morphology  = 'in.swc'
metrics = 'out.json'

m = Morpho(morphology)

print(title)
print('neuron:', morphology)

all_nodes = [i for i in m.traverse()]
all_branches = [b for b in m.branches()]
morphometrics = {}

neurites = set([m.neurite(i) for i in all_nodes])
print('neurites:', ' '.join(n for n in neurites))

for neurite in neurites:
    measures = {}
    branches = filter(lambda b: m.neurite(b[0])==neurite, all_branches)
    nodes = filter(lambda i: m.neurite(i)==neurite, all_nodes)
    bifurcations = filter(lambda i: m.is_fork(i), nodes)
    tips = filter(lambda i: m.is_leaf(i), nodes)
    stems = filter(lambda b: m.order(b[0])==1, branches)

    area = sum(m.area(i) for i in nodes)
    length = sum(m.length(i) for i in nodes)
    volume = sum(m.volume(i) for i in nodes)
    measures['area'] = area
    measures['length'] = length
    measures['volume'] = volume
    diams =[m.diam(i) for i in nodes]
    measures['extent euclidean'] = [s for s in m.size(idents=nodes)]
    measures['extent radial'] = max([m.distance(i, radial=True) for i in nodes])
    measures['extent path'] = max([m.distance(i) for i in tips])
    measures['diameters'] = np.mean(diams), np.min(diams), np.max(diams)
    measures['diameter effective'] = area/(pi*length)
    measures['order'] = max([m.order(i) for i in tips])
    measures['stems'] = len(stems)
    measures['bifurcations'] = len(bifurcations)
    measures['tips'] = len(tips)
    measures['branches'] = len(branches)
    morphometrics[neurite] = measures

with open(metrics, 'w') as file:
    json.dump(morphometrics, file, indent=4, sort_keys=True)
print('morphometric data saved to', metrics)
