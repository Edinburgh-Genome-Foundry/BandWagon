from .BandsPattern import Band, BandsPattern


def generate_ladder(bands_migrations):
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
    ])

LADDER_100_to_4k = generate_ladder({
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
