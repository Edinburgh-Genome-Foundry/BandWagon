GelSimulator.py
================

Gel Simulator is a Python library to predict gel migration patterns
from enzymatic DNA digestions. It supports hundreds of enzymes (thanks to BioPython),
single- and multiple-enzymes digestions, and custom ladders.

It uses Matplotlib to produce plots like this one:

.. figure:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GelSimulator.py/master/examples/multi_enzyme_mixes.png
    :align: center

While fairly minimal, Gel Simulator was written so as to be easily reusable, customizable, and extensible.

License = MIT
---------------

Gel Simulator is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/GelSimulator.py>`_ under the MIT licence.

Everyone is welcome to contribute !

Installation
--------------

With PIP:

.. code:: python

    (sudo) pip install gelsimulator

Gel Simulator can be installed by unzipping the source code in one directory and using this command:

.. code:: python

    (sudo) python setup.py install

PIP install is coming soon !



Examples of use
----------------


Computing digestion bands sizes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This first example shows how to compute digestion bands in the case of
a linear fragment, a circular fragment, and a multi-enzymes digestion.

..  code:: python

    from gelsimulator import compute_digestion_bands

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

Plotting a gel simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~

The following example shows how to plot the the digestion patterns produced
by different restriction enzymes on a same DNA sequence:

..  code:: python

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

.. figure:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/GelSimulator.py/master/examples/simple_gel_simulation.png
    :align: center

Using a custom ladder
~~~~~~~~~~~~~~~~~~~~~~~

In the previous example, we use the pre-defined `LADDER_100_to_4k`, but what
if we wanted another ladder that Gel Simulator does not provide ?

New ladders can be easily defined by providing a dict of the form:

..  code:: python

    {
    `band_size_1: y_coordinate_1`,
    `band_size_2: y_coordinate_2`,
    `band_size_3: y_coordinate_3`
    }

Where `band_size` is the known length of the ladder's DNA fragment,
and `y_coordinate` the y coordinate in pixels of the corresponding band in a
picture of the ladder's migration.

Note that you must provide at least 3-4 bands for the ladder to be meaningful.
The gel must be oriented with larger-fragments bands on top (as is usually the case)
It is not necessary to provide any "origin" of the ladder as it will be
computed automatically. For instance:

..  code:: python

    custom_ladder = GelLadder(bands={
        # band_size : meaured y-coordinate
        100: 200,
        300: 170,
        500: 150,
        1650: 100,
        4000: 65
    })
    gel_simulator = GelSimulator(custom_ladder)


See file `examples/gel_simulation_with_enzymes_mixes.py` for a more complete
example involving a custom ladder and multi-enzyme digestions.
