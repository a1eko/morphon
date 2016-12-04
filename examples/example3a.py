from __future__ import print_function, division
from morphon import Morpho

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


title = 'Neuron 3D plot'
morphology  = 'in2.swc'

print(title)
print('neuron:', morphology)

m = Morpho(morphology)
fig = plt.figure()
ax = fig.gca(projection='3d')
for branch in m.branches():
    x = [m.coord(i)[0] for i in branch]
    y = [m.coord(i)[1] for i in branch]
    z = [m.coord(i)[2] for i in branch]
    parent = m.parent(branch[0])
    if parent is not None:
        x.insert(0, m.coord(parent)[0])
        y.insert(0, m.coord(parent)[1])
        z.insert(0, m.coord(parent)[2])
    ax.plot(x, z, y, color='black')

plt.show()
