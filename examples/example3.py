from __future__ import print_function, division
from morphon import Morpho, plot

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


title = 'Neuron 3D plot'
morphology  = 'in.swc'

print(title)
print('neuron:', morphology)

m = Morpho(morphology)
fig = plt.figure()
ax = fig.gca(projection='3d')
plot(m, ax, projection='3d')

plt.show()
