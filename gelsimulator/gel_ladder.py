import numpy as np

class GelLadder:
    """Class to represent gel ladders. These ladders serve as a scale for
    plotting any other gel simulation."""

    def __init__(self, bands):
        self.bands = bands
        self.compute_migration_distance_parameters()

    def compute_migration_distance_parameters(self):
        """Fit the ladder's band with a mathematical model

        This allows to then predict how far new bands will migrate with respect
        to this ladder.
        """
        sizes, distances = [
            np.array(a)
            for a in
            zip(*sorted(self.bands.items()))
        ]

        # This very hackish function fitting gets rid of the scipy dependency
        (a_est, b_est), r, _, _, _ = min(
            [
                np.polyfit(sizes, np.log(distances-c), 1, full=True)
                for c in np.linspace(0, 0.9*min(distances), 10)
            ],
            key = lambda abr: (abr[1]**2).sum()
        )
        self.parameters = a_est, b_est


    def band_size_to_migration(self, band_size):
        """Compute how far a fragment of size `band_size` will migrate.

        This uses a mathematical fit of the ladder's bands.
        """
        a, b = self.parameters
        return np.exp(a * band_size + b)

    def migration_to_band_size(self, migration):
        a, b = self.parameters
        return (np.log(migration) - b) / a


    def migration_distances_span(self):
        min_band = min(self.bands.keys())
        max_band = max(self.bands.keys())
        return [self.band_size_to_migration(band)
                for band in (max_band, min_band)]


LADDER_100_to_4k = GelLadder(bands={
    100: 205,
    200: 186,
    300: 171,
    400: 158,
    500: 149,
    650: 139,
    850: 128,
    1000: 121,
    1650: 100,
    2000: 90,
    3000: 73,
    4000: 65
})
