from __future__ import print_function, division
from morphon import Morpho, measure
import json


title = 'Basic morphometrics'
morphology  = 'in.swc'
metrics = 'out.json'

m = Morpho(morphology)

print(title)
print('neuron:', morphology)

morphometrics = measure(m)

with open(metrics, 'w') as file:
    json.dump(morphometrics, file, indent=4, sort_keys=True)
print('morphometric data saved to', metrics)
