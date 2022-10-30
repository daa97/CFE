import numpy as np
from scipy.interpolate import CubicSpline
import csv

rang = np.array([[1, 1.1, 1.2, 1.32, 1.48, 1.65, 1.8, 2, 2.25, 2.5, 2.8, 3.1, 3.5, 3.9, 4.4, 4.9, 5.4, 6, 6.7, 7.4, 8.25, 9.1]])
P = rang * np.array([[10**x for x in range(2,8)]]).transpose()
P = np.round(P, 6)
P = P.flatten()
T = np.arange(20, 205, 10, dtype=int)



table_in = "CEA Property Tables\\s cryo.csv"

with open(table_in, mode='r') as f:
    r = csv.reader(f)
    contents = np.array([row for row in r])
header = contents[0,1:].astype(float) * 101325
left = contents[1:,0].astype(int)

data = contents[1:,1:].astype(float)
print(data)
splines = []
#splinetemps = []
for i in range(len(left)):
    row = data[i,:]
    di = row[np.isfinite(row)]
    hi = header[np.isfinite(row)]
    try:
        splines.append(CubicSpline(hi, di,extrapolate=False))
    except:
        splines.append(np.nan)


values = np.zeros((len(T),len(P)))
for i in range(len(T)):
    a = np.where(left==T[i])[0]
    print(a, T[i])
    try:
        s = splines[a[0]]
    except:
        pass
    for j in range(len(P)):
        
        try:
            values[i,j] = s(P[j])
            print(s(P[j]))
            print("YAY", values[i,j])
        except:
            print("ouch")
            values[i,j] = np.nan
    print(values)
print(values.tolist())
with open("out.csv", mode='w', newline='') as f:
    a = csv.writer(f)
    a.writerow(P)
    a.writerows(values)
    