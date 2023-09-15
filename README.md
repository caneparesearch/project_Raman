# project_Raman
Workflow code and output for computing Raman spectroscopy and properties related to lattice vibration using CRYSTAL.
## cryoutput_parsing_databse
*crystal17_output_files* contains all CRYSTAL output files (.out) that include all the calculated thermodynamic and vibrational properties.\
*crystal_out.py* parses the CRYSTAL output files.\
*crystal_mongo.py* writes to and reads from our MongoDB database.
## web_interface
Code to build a web interface that allows users to interactively query compounds and visualize their Raman spectra.
