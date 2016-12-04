from __future__ import print_function, division
from morphon import Morpho


title = 'Created morphology'
morphology  = 'out.swc'

m = Morpho()
root = m.add('soma', [0, 0, 0], 10)

node = m.add('dend', [10, 0, 0], 2, parent=root)
node = m.add('dend', [20, 0, 0], 1, parent=node)
fork = m.add('dend', [30, 0, 0], 0.8, parent=node)
node = m.add('dend', [35, 5, 0], 0.7, parent=fork)
node = m.add('dend', [40, 10, 0], 0.7, parent=node)

node = m.add('axon', [0, 10, 0], 1, parent=root)
node = m.add('axon', [0, 50, 0], 0.5, parent=node)

print(title)
print(m)
m.save(morphology, header=title)
print('saved to', morphology)
