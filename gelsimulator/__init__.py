""" gelsimulator/__init__.py """

# from .GelSimulator import GelSimulator
# from .GelLadder import GelLadder, LADDER_100_to_4k
from ladders import LADDER_100_to_4k
from BandsPattern import Band, BandsPattern, BandsPatternsSet
from .tools import (merge_bands_in_pattern,
                    bands_are_similar,
                    bands_patterns_are_similar,
                    compute_digestion_bands)
from .version import __version__
