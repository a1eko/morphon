from __future__ import print_function, division
from morphon import Morpho

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


title = 'Neuron 3D plot'
morphology  = 'in.swc'

print(title)
print('neuron:', morphology)

m = Morpho(morphology)
fig = plt.figure()
ax = fig.gca(projection='3d', title=morphology)

all_nodes = [i for i in m.traverse()]
all_branches = [b for b in m.branches()]
neurites = set([m.neurite(i) for i in all_nodes])
branches = {}
for neurite in neurites:
    branches[neurite] = filter(lambda b: m.neurite(b[0])==neurite, all_branches)

def plot(branches, ax, color='black', linewidth=1):
    for branch in branches:
        x = [m.coord(i)[0] for i in branch]
        y = [m.coord(i)[1] for i in branch]
        z = [m.coord(i)[2] for i in branch]
        parent = m.parent(branch[0])
        if parent is not None:
            x.insert(0, m.coord(parent)[0])
            y.insert(0, m.coord(parent)[1])
            z.insert(0, m.coord(parent)[2])
        ax.plot(x, z, y, color=color, linewidth=linewidth)

color = {'soma': 'black', 'dend': 'blue', 'apic': 'green', 'axon': 'red'}
thickness = {'soma': 5, 'dend': 1, 'apic': 1, 'axon': 0.5}

for neurite in neurites:
    plot(branches[neurite], ax, color=color[neurite], linewidth=thickness[neurite])

plt.show()
