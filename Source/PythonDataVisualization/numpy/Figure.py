
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-3, 3, 50)

y1 = 2*x +1
y2 = x**2


plt.figure(num = 1, figsize=(8, 5))
plt.plot(x, y2)
plt.plot(x, y1, color='red', linewidth=1.0, linestyle='--')
plt.show()



