"""Simple gel simulation with gelsimulation.py.

This example shows how to plot the the digestion patterns produced by different
restriction enzymes on a same DNA sequence.
"""

from gelsimulator import GelSimulator, LADDER_100_to_4k
import matplotlib.pyplot as plt

with open("example_sequence.txt", "r") as f:
    sequence = f.read()

enzymes = ["BamHI", "EcoRI", "EcoRV", "PstI", "SpeI", "XbaI"]
gel_simulator = GelSimulator(LADDER_100_to_4k)

fig, ax = plt.subplots(1, figsize=(1.1*(len(enzymes)+1), 5))
gel_simulator.format_ax(ax)
gel_simulator.plot_ladder(ax, x_coord=1)
for i, enzyme in enumerate(enzymes):
    gel_simulator.plot_digestion_result(ax, sequence, [enzyme],
                                        x_coord=i+2, label=enzyme)

fig.savefig("gel_simulation.png", bbox_inches="tight")
