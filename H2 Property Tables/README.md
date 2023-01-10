# Hydrogen Property Tables
These scripts are used to create csv tables containing property data for hydrogen across a range of temperatures. The script `refprop.py` is used to collect data from NIST Refprop up to 1500 K. The script `CEA.ipynb` is used to collect data from NASA CEA above 1000 K. The script `combine_data.ipynb` combines these datasets into a single numpy archive `H2_compiled.npz` and exports them as a set of CSVs.

The tables which are copied into the `/Updated` folder will then be referenced by `fluids.py` in the repository root folder.
