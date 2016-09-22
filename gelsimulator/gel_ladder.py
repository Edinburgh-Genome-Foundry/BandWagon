import numpy as np
from scipy.optimize import curve_fit


class GelLadder:
    """Class to represent gel ladders. These ladders serve as a scale for
    plotting any other gel simulation."""

    def __init__(self, bands):
        self.bands = bands
        self.compute_migration_distance_predictor()

    def compute_migration_distance_predictor(self):

        sizes, distances = [
            np.array(a)
            for a in
            zip(*sorted(self.bands.items()))
        ]
        a_est, b_est = np.polyfit(sizes, np.log(distances), 1)

        def predict_migration_distance(size, a, b, c):
            return np.exp(a * size + b) + c
        (a, b, c), _ = curve_fit(
            predict_migration_distance,
            sizes,
            distances,
            [a_est, b_est, 0]
        )

        def migration_predictor(size):
            return predict_migration_distance(size, a, b, c) - c
        self._migration_predictor = migration_predictor

    def compute_migration_distance(self, band_size):
        return self._migration_predictor(band_size)

    @property
    def migration_distances_span(self):
        min_band = min(self.bands.values())
        max_band = max(self.bands.values())
        return [self.compute_migration_distance(band)
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
