from bandwagon import (plot_all_digestion_patterns, plot_records_digestions,
                       load_record, LADDER_100_to_4k)
import os

records = [
    load_record(os.path.join('records', filename),
                id=filename, topology='circular')
    for filename in sorted(os.listdir('records'))
]
digestions= [
    ('BamHI', 'NcoI'),
    ('BsaI', 'XbaI'),
    ('StyI',)
]

plot_records_digestions(
    records=records,
    digestions=digestions,
    ladder=LADDER_100_to_4k,
    target="plot_records_digestions_example.pdf")

axes = plot_all_digestion_patterns(
    records=records,
    digestions=digestions,
    ladder=LADDER_100_to_4k)
axes[0].figure.savefig("plot_all_digestion_patterns.png", bbox_inches='tight')