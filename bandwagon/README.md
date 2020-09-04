# Code organization

This document walks you through the BandWagon code. Please request changes if anything is unclear.

**Core classes**

- ``Band`` (defined in ``Band.py``) is a class representing a single band,
  with its color, size, label, etc.
- ``BandsPattern`` (defined in ``BandsPattern.py``) is a class representing a
  pattern comprising several bands. The pattern may also have a top label,
  a corner note, an image (for instance, of a gel) on the side, etc.
- ``BandsPatternsSet`` (defined in ``BandsPatternsSet.py``) is a class
  representing a list of BandPatterns which will be plotted side by side in
  a row, with possibly a vertical label on the left side, a ladder, etc.

**Other modules**

- ``ladders.py`` implements some ways to define a ladder (which is really just a dict) either on the user side, or from an AATI Fragment Analyser calibration table.
-  ``plot_records_digestions.py`` implements methods allowing to plot a bands pattern along with a schema of the sequence indicating where the digestion cuts.
-  ``tools.py`` implements generic methods for handling biopython records, etc.
