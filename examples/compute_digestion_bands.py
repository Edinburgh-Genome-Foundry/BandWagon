"""Prediction of digestion bands sizes in different scenarios.

The following example shows how to compute digestion bands in the case of
a linear fragment, a circular fragment, and a multi-enzymes digestion.
"""

from bandwagon import compute_digestion_bands

# Read the sequence (a string of the form 'ATGTGTGGTA...' etc.)
with open("example_sequence.txt", "r") as f:
    sequence = f.read()

# Compute digestion bands for a linear construct
print(compute_digestion_bands(sequence, ["EcoRI"], linear=True))
# Result >>> [400, 1017, 3583]

# Compute digestion bands for a circular construct
print(compute_digestion_bands(sequence, ["EcoRI"], linear=False))
# Result >>> [1017, 3983]

# Compute digestion bands for an enzymatic mix
print(compute_digestion_bands(sequence, ["EcoRI", "BamHI"]))
# Result >>> [400, 417, 600, 3583]
