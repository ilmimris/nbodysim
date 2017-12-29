import matplotlib.pyplot as plt
import numpy as np
import time
start_time = time.time()

n = 3
dt = 0.01
t_end = 10
m = np.zeros(n)
m[:] = 1

r = np.zeros([n,3])
v = np.zeros([n,3])
a = np.zeros([n,3])
jk = np.zeros([n,3])

def elips():
	for i in range(n):
		phi = i * 2 * np.pi / n
		r[i,0] = np.cos(phi)
		r[i,1] = np.sin(phi)
		r[i,2] = 0
	v_abs = 1.0/np.sqrt(np.sqrt(3))
	for i in range(n):
		phi = i * 2 * np.pi / n
		v[i,0] = - v_abs * np.sin(phi)
		v[i,1] = v_abs * np.cos(phi)
		v[i,2] = 0

def circle():
	for i in range(n):
		phi = i * 2 * np.pi / n
		r[i,0] = np.cos(phi)
		r[i,1] = np.sin(phi)
		r[i,2] = 0
	v_abs = .5/np.sqrt(np.sqrt(3))
	for i in range(n):
		phi = i * 2 * np.pi / n
		v[i,0] = - v_abs * np.sin(phi)
		v[i,1] = v_abs * np.cos(phi)
		v[i,2] = 0

def inf():
	## Inf shape
	r[0][0] = 0.9700436
	r[0][1] = -0.24308753
	r[0][2] = 0
	v[0][0] = 0.466203685
	v[0][1] = 0.43236573
	v[0][2] = 0

	r[1][0] = -r[0][0]
	r[1][1] = -r[0][1]
	r[1][2] = -r[0][2]
	v[1][0] = v[0][0]
	v[1][1] = v[0][1]
	v[1][2] = v[0][2]

	r[2][0] = 0
	r[2][1] = 0
	r[2][2] = 0
	v[2][0] = -2 * v[0][0]
	v[2][1] = -2 * v[0][1]
	v[2][2] = -2 * v[0][2]
	## End Shape inf

def nol_V():
	for i in range(n):
		phi = i * 2 * np.pi / n
		r[i,0] = np.cos(phi)
		r[i,1] = np.sin(phi)
		r[i,2] = 0
	v_abs = 0.
	for i in range(n):
		phi = i * 2 * np.pi / n
		v[i,0] = - v_abs * np.sin(phi)
		v[i,1] = v_abs * np.cos(phi)
		v[i,2] = 0

elips()

for i in range(n):
	for k in range(3):
		a[i,k] = jk[i,k] = 0.0

for i in range(n):
	j = i+1
	for j in range(j,n):
		rji = np.zeros(3)
		vji = np.zeros(3)

		for k in range(3):
			rji[k] = r[j,k] - r[i,k]
			vji[k] = v[j,k] - v[i,k]
			
		r2 = 0
		for k in range(3):
			r2 += rji[k] * rji[k]
		
		r3 = r2 * np.sqrt(r2)
		
		rv = 0
		for k in range(3):
			rv += rji[k] * vji[k]
		
		rv /= r2

		for k in range(3):
			a[i,k] += m[j] * rji[k] / r3
			a[j,k] -= m[i] * rji[k] / r3
			jk[i,k] += m[j] * (vji[k] - 3 * rv * rji[k]) / r3
			jk[j,k] -= m[i] * (vji[k] - 3 * rv * rji[k]) / r3			
		
	


v[0,0] += 0.0001
ekin = 0
epot = 0

for i in range(n):
	j = i + 1
	
	for j in range(j,n):
		rji = np.zeros(3)

		for k in range(3):
			rji[k] = r[j,k] - r[i,k]

		r2 = 0
		for k in range(3):
			r2 += rji[k] * rji[k]

		r_ = np.sqrt(r2)
		epot -= m[i]*m[j]/r_

	for k in range(3):
		ekin += 0.5 * m[i] * v[i,k] * v[i,k]

e_in = ekin + epot
print ("Initial total energy E_in = {} \n".format(e_in))

dt_out = 0.01
t_out = dt_out

old_r=np.zeros([n,3])
old_v=np.zeros([n,3])
old_a=np.zeros([n,3])
old_j=np.zeros([n,3])

for t in np.arange(0,t_end,dt):
	for i in range(n):
		for k in range(3):
			old_r[i,k] = r[i,k]
			old_v[i,k] = v[i,k]
			old_a[i,k] = a[i,k]
			old_j[i,k] = jk[i,k]
			r[i,k] += v[i,k]*dt + a[i,k]*dt*dt/2 + jk[i,k]*dt*dt*dt/6
			v[i,k] += a[i,k]*dt + jk[i,k]*dt*dt/2

	for i in range(n):
		for k in range(3):
			a[i,k] = jk[i,k] = 0.0

	for i in range(n):
		j = i + 1
		for j in range(j,n):
			rji = np.zeros(3)
			vji = np.zeros(3)
		
			for k in range(3):
				rji[k] = r[j,k] - r[i,k]
				vji[k] = v[j,k] - v[i,k]
			r2 = 0
			for k in range(3):
				r2 += rji[k] * rji[k]
			r3 = r2 * np.sqrt(r2)
			rv = 0			
			for k in range(3):
				rv += rji[k] * vji[k]
			rv /= r2
			
			for k in range(3):
				a[i,k] += m[i] * rji[k] / r3
				a[j,k] -= m[i] * rji[k] / r3
				jk[i,k] += m[i] * (vji[k] - 3 * rv * rji[k]) / r3
				jk[j,k] -= m[i] * (vji[k] - 3 * rv * rji[k]) / r3

	for i in range(n):
		for k in range(3):
			r[i,k] = old_r[i,k] + (old_v[i,k] + v[i,k])*dt/2 + (old_a[i,k] - a[i,k])*dt*dt/12
			v[i,k] = old_v[i,k] + (old_a[i,k] + a[i,k])*dt/2 + (old_j[i,k] - jk[i,k])*dt*dt/12
			

	# if (t >= t_out):
	# 	for i in range(n):
	# 		for k in range(3):
	# 			print  " {} ".format(r[i,k])
	# 		for k in range(3):
	# 			print " {} ".format(v[i][k])
	# 		print "\n"
	# 	t_out += dt_out
	plt.plot(r[:,0],r[:,1],"r.")

epot = ekin = 0
for i in range(n):
	j = i + 1
	for j in range(j,n):
		rji=np.zeros(3)

		for k in range(3):
			rji[k] = r[j,k] - r[i,k]
		r2 = 0
		
		for k in range(3):
			r2 += rji[k] * rji[k]

		epot -= m[i]*m[i]/np.sqrt(r2)

	for k in range(3):
		ekin += 0.5 * m[i] * v[i,k] * v[i,k]


e_out = ekin + epot
print ("Final total energy E_out = {}".format(e_out))
print ("absolute energy error: E_out - E_in = {}".format(e_out - e_in))
print ("relative energy error: (E_out - E_in) / E_in = {}".format((e_out - e_in) / e_in ))
print("Finished %s seconds" % (time.time() - start_time))
plt.xlabel("X (AU)")
plt.ylabel("Y (AU)")
plt.show()
