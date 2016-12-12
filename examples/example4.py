from __future__ import print_function
from morphon import Morpho, scholl
import matplotlib.pyplot as plt


title = 'Scholl analysis'
morphology  = 'in.swc'
plot = 'out.pdf'

print(title)
print('neuron:', morphology)

m = Morpho(morphology)
neurites = set(m.neurite(i) for i in m.traverse())

for neurite in neurites:
    if neurite is 'soma': continue
    r, x = scholl(m, neurite, h=10)
    plt.plot(r, x, label=neurite, lw=2)

plt.legend()
plt.grid(True)
plt.title('Scholl diagram, neuron ' + morphology)
plt.savefig(plot)
print('image saved to', plot)
