# YANKsuite

from YankSuite import *

saves data and images to the output directory of the yaml script passed in.

if sim1, sim2, and sim3 are each repeats of the same YANK experiment\
sim1 = YankAnalyzer("full/path/to/yaml/that/ran/YANK/simulation1.yaml")\
sim2 = YankAnalyzer("full/path/to/yaml/that/ran/YANK/simulation2.yaml")\
sim3 = YankAnalyzer("full/path/to/yaml/that/ran/YANK/simulation3.yaml")\
sims = YankMultiAnalyzer({'data':[sim1, sim2, sim3]})

try
for indv simulations
sim1.extract_meaningful_trajs() # Extracts physically meaningful trajs/states 