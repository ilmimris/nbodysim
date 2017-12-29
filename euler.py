"""
30-Oktober-2017 Prodi Astronomi
2-Body simulation with euler method

"""

import matplotlib.pyplot as plt 
import numpy as np
import time
start_time = time.time()

r=[1,0,0]
v=[0,0.5,0]
a=np.zeros(3)

# dt=float(input("Please provide a value for the time step : ")) ## masukkan 0.001
dt = 0.0001

ekin = 0.5 * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
epot = -1.0/np.sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])
e_in = ekin + epot

print ("Intial total energi E_in = {}".format(e_in))

dt_out = 0.01
t_out = dt_out

# t_end = int(input("time stop at = "))  ## masukkan 10
t_end = 10

for t in np.arange(0,t_end,dt):
	r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]
	for k in range(3):
		a[k] = - r[k] / (r2 * np.sqrt(r2))
	for k in range(3):
		r[k] += v[k] * dt
		v[k] += a[k] * dt
	# if t >= t_out:
	# 	print "{} {} {} ".format(r[0],r[1],r[2])
	# 	print "{} {} {} \n".format(v[0],v[1],v[2])
	# 	t_out += dt_out
	plt.plot(r[0],r[1],"r.")		

ekin = 0.5 * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
epot = -1.0/np.sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])
e_out = ekin + epot
print ("Final total energy E_out = {} \n".format(e_out))
print ("absolute energy error: E_out - E_in = {} \n".format(e_out - e_in))
print ("relative energy error: (E_out - E_in) / E_in = {} \n".format((e_out - e_in) / e_in))
print("Finished in %s seconds" % (time.time() - start_time))
plt.xlabel("X (AU)")
plt.ylabel("Y (AU)")
plt.show()