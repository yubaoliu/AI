import numpy as np



open, close, volume =np.loadtxt('000001.csv', delimiter=',', skiprows=1, usecols=(1, 4, 6), unpack=True)


change = close-open


print(change)
