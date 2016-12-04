from __future__ import print_function, division
from morphon import Morpho
import matplotlib.pyplot as plt
import numpy as np


title = 'Scholl analysis'
morphology  = 'in2.swc'
plot = 'out4.pdf'
h = 10.0

def crossings(m, i, h):
    p = m.parent(i)
    r1 = m.distance(p, radial=True)
    r2 = m.distance(i, radial=True)
    k0 = int(r1/h)
    kn = int(r2/h)
    return kn-k0, k0

def scholl(m, neurite, h):
    nodes = filter(lambda i: m.neurite(i)==neurite, all_nodes)
    rmax = max([m.distance(i, radial=True) for i in nodes])
    radx = np.array([k*h for k in range(int(rmax/h))])
    crox = np.zeros(int(rmax/h), dtype=int)
    for i in nodes:
        ncross, icross = crossings(m, i, h)
        if ncross > 0:
            for k in range(ncross):
                crox[icross+k] += 1
        elif ncross < 0:
            for k in range(-ncross):
                crox[icross-k] += 1
    return radx, crox


print(title)
print('neuron:', morphology)

m = Morpho(morphology)
all_nodes = [i for i in m.traverse()]
neurites = set([m.neurite(i) for i in all_nodes])

for neurite in neurites:
    if neurite is 'soma': continue
    r, x = scholl(m, neurite, h)
    plt.plot(r, x, label=neurite, lw=2)

plt.legend()
plt.grid(True)
plt.title('Scholl diagram, neuron ' + morphology)
plt.savefig(plot)
print('image saved to', plot)
