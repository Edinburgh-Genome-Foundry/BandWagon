
class GelLadder:

    def __init__(self, bands, use_log_interpolation=True):
        self.bands = bands
        sizes, distances = [
            np.array(a)
            for a in
            zip(*sorted(bands.items()))
        ]
        a_est, b_est = np.polyfit(sizes, np.log(distances), 1)

        def size_to_migration_distance(s, a, b, c):
            return np.exp(a * s + b) + c
        (a, b, c), _ = curve_fit(
            size_to_migration_distance,
            sizes,
            distances,
            [a_est, b_est, 0]
        )

        def migration_predictor(s):
            return size_to_migration_distance(s, a, b, c) - c
        self._migration_predictor = migration_predictor

    def compute_migration_distance(self, band_size):
        return self._migration_predictor(band_size)

    @property
    def migration_distances_span(self):
        min_band = min(self.bands.values())
        max_band = max(self.bands.values())
        return [self.compute_migration_distance(band)
                for band in (max_band, min_band)]
