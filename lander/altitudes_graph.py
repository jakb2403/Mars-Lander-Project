import numpy as np
import matplotlib.pyplot as plt
results = np.loadtxt('altitudes.txt')
K_h = 0.025
target_descent_rate = -(0.5 + K_h * (results[:,0]))
print(results.shape)
print(target_descent_rate.shape)
plt.figure(1)
plt.clf()
plt.ylabel('Descent rate (m/s)')
plt.xlabel('Altitude (m)')
plt.grid()
plt.plot(results[:, 0], results[:, 1], label='Actual descent rate')
plt.plot(results[:, 0], target_descent_rate, label='Target descent rate')
plt.legend()
plt.show()