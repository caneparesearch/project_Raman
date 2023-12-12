# project_Raman
Workflow code and output for computing Raman spectroscopy and properties related to lattice vibration using CRYSTAL.
## cryoutput_parsing_databse
*crystal_output_files* contains all CRYSTAL output files (.out) that include all the calculated thermodynamic and vibrational properties.\
*crystal_out.py* parses the CRYSTAL output files.\
*crystal_mongo.py* writes to and reads from our MongoDB database.
## web_interface
Code to build a web app (https://raman-db.streamlit.app/) that allows users to search for chemical formula, select the desired compound according to its ICSD id, view the crystal
structure, and interactively plot the Raman and IR spectra with different convolution schemes.

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
