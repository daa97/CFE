import numpy as np, CEAio

P = np.array([[round(10**x,5) for x in range(2,9,1/30)]]).transpose()
T = np.arange(500, 8e3, 5, dtype=int)

files = CEAio.write_inp(P,T)
CEAio.run_cea(files)
props = CEAio.read_out(files)
props['VOLUME'] = 1/np.array(props["RHO, KG/CU M"])
tables = CEAio.tabulate(props, P, T)
np.savez("CEAprops.npz", tables)

#CEAio.write_csvs(tables)