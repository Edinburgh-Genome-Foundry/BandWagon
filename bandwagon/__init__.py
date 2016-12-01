""" bandwagon/__init__.py """

# from .BandWagon import BandWagon
# from .GelLadder import GelLadder, LADDER_100_to_4k
from .ladders import LADDER_100_to_4k, generate_ladder
from .BandsPattern import Band, BandsPattern, BandsPatternsSet
from .tools import (merge_bands_in_pattern,
                    bands_are_similar,
                    bands_patterns_are_similar,
                    compute_digestion_bands)
from .version import __version__
