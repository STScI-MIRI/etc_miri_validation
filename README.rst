ETC Validation for MIRI Imager
==============================

This repository contains code to run the JWST ETC (via the Pandeia engine) and compare the predictions for signal-to-noise ratio and background surface brightness to observations.

This project was finished in November 2024.


Description
-----------
The goal of this project was to validate the accuracy of the JWST Exposure Time Calculator (ETC) (version 3.2) for the MIRI Imager. This was done by comparing the predicted signal-to-noise ratio (S/N) and background surface brightness from the ETC to the values measured from the absolute flux calibration data for three A dwarfs, using aperture photometry measurements. In most cases, the ETC 3.2 (which does not include corrections for imager count-rate loss or thermal background changes) overestimates the S/N and underestimates the background level, but the discrepancies depend on the selected star and filter.

The methods and results of this project are described in detail in SOCCER document JWST-STScI-008884.


Contributors
------------

Marjorie Decleir


License
-------

This project is Copyright (c) Marjorie Decleir and licensed under
the terms of the BSD 3-Clause License (see the ``LICENSE`` file for more information).


Dependencies
------------

This project uses the `JWST ETC Pandeia Engine <https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-pandeia-engine-tutorial#gsc.tab=0>`_ to obtain ETC predictions for signal-to-noise ratio and background surface brightness.


Citation
--------
If this repository was useful for your work, please cite the SOCCER document JWST-STScI-008884.
