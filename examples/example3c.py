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

neurites = set([m.neurite(i) for i in m.traverse()])

def plot(sections, ax, color=None, linewidth=None):
    for section in sections:
        x, y, z = m.coords(section)
        ax.plot(x, z, y, color=color, linewidth=linewidth)

color = {'soma': 'black', 'dend': 'blue', 'apic': 'green', 'axon': 'red'}
thickness = {'soma': 5, 'dend': 1, 'apic': 1, 'axon': 0.5}

for neurite in neurites:
    plot(m.sections(with_parent=True, neurites=[neurite]), 
        ax, color=color[neurite], linewidth=thickness[neurite])

plt.show()
