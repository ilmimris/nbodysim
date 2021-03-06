import matplotlib.pyplot as plt
import numpy as np
import time
start_time = time.time()

r=[1,0,0]
v=[0,0.5,0]
a=np.zeros(3)

# dt=float(input("Please provide a value for the time step : ")) ## masukkan 0.001
dt = 0.001

ekin = 0.5 * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
epot = -1.0/np.sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])
e_in = ekin + epot

print ("Intial total energi E_in = {}".format(e_in))

r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
for k in range(3):
	a[k] = - r[k] / (r2 * np.sqrt(r2));

dt_out = 0.01;
t_out = dt_out;

# t_end = int(input("time stop at = "))  ## masukkan 10
t_end = 10

for t in np.arange(0,t_end,dt):
	for k in range(3):
		v[k] += 0.5 * a[k] * dt
		r[k] += v[k] * dt
	r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]
	for k in range(3):
		a[k] = - r[k] / (r2 * np.sqrt(r2))
		v[k] += 0.5 * a[k] * dt
	# if (t >= t_out):
	# 	print "{} {} {} ".format(r[0],r[1],r[2])
	# 	print "{} {} {} \n".format(v[0],v[1],v[2])
	# 	t_out += dt_out;
	plt.plot(r[0],r[1],"r.",linewidth=0.1)

ekin = 0.5 * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
epot = -1.0/np.sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
e_out = ekin + epot;

print ("Final total energy E_out = {}".format(e_out))
print ("absolute energy error: E_out - E_in = {}".format(e_out - e_in))
print ("relative energy error: (E_out - E_in) / E_in = {}".format((e_out - e_in) / e_in))
print("Finished in %s seconds" % (time.time() - start_time))
plt.xlabel("X (AU)")
plt.ylabel("Y (AU)")
plt.show()