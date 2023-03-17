import numpy as np
import matplotlib.pyplot as plt

betas = np.linspace(-np.pi/2,np.pi/2,100)
L = np.zeros((len(betas),))
for i,beta_4_rad in enumerate(betas):
    if beta_4_rad + np.pi/6 > 0:
        n = 3 
    else:
        n = 2
    L[i] = 0.5 * (np.abs(np.sin(beta_4_rad + np.pi/6)))**n
betas = np.multiply(betas,180/np.pi)
plt.plot(betas,L)
plt.show()