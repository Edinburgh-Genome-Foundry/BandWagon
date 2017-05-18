"""Definition of custom ladders for BandWagon."""

from .BandsPattern import Band, BandsPattern
try:
    import pandas
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

def custom_ladder(label, bands_migrations):
    """Generate a BandsPattern object meant to be used as a ladder.

    By default the ladder is dark red on a white background.
    """
    migrations = bands_migrations.values()
    min_migrations = min(migrations)
    ladder_span = max(migrations) - min_migrations
    if min_migrations < 0.1*ladder_span:
        shift = 0
    else:
        shift = min_migrations - 0.1*ladder_span
    return BandsPattern([
        Band(size, migration_distance=migration - shift, band_color="#8B0000")
        for (size, migration) in bands_migrations.items()
    ], background_color="#ffffff", label=label)


def ladder_from_aati_fa_calibration_table(filepath=None, dataframe=None,
                                          label=None):
    """Extract a BandPattern from an AATI Fragment Analyzer calibration file.

    The calibration table is generated after each run from the migration
    times of the ladder and is generally in a file called something like
    ``2016 11 08 16H 57M Size Calibration.csv``.

    This method requires Pandas installed.

    Parameters
    ----------

    filepath
      Path to the calibration file.

    dataframe
      Pandas dataframe obtained by reading the file, can be provided instead of
      ``filepath``

    label
      Label that will be given to the ladder when plotted.
    """
    if not PANDAS_AVAILABLE:
        raise ImportError("The Pandas library is required for this method.")

    if filepath is not None:
        dataframe = pandas.read_csv(filepath)

    dataframe["migration"] = (1.1 * dataframe["Time (sec)"].max() -
                              dataframe["Time (sec)"])

    return custom_ladder(label, {
        row["Ladder Size (bp)"]: row["migration"]
        for i, row in dataframe.iterrows()
    })

LADDER_100_to_4k = custom_ladder("100-4k", {
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
