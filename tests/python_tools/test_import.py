import sys
import numpy as np
sys.path.append('/storage/gpfs_data/neutrino/users/gi/sand-reco/tests/python_tools')
from NLLtreeoutReader import Reader
from Helix import Helix
from EventDisplay import EventDisplay
import matplotlib.pyplot as plt
test_file = "/storage/gpfs_data/neutrino/users/gi/sand-reco/tests/test_reconstruct_NLLmethod_smear.root"
reader = Reader(test_file, "tReco")
my_helix_true = reader.get_true_helix(9)
wires = reader.get_wire_info(9)

fig, ax = plt.subplots(1, 2, figsize=(20, 10))

display = EventDisplay(arg_file_name = test_file, 
                       arg_event_idx = 9, 
                       arg_true_helix_points = my_helix_true.get_helix_points(), 
                       ax_0 = ax[0], 
                       ax_1 = ax[1],
                       arg_wires_info = wires)

display.plot_sand()