import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
start_time = time.time()

n = 3
dt = 0.01
t_end = 1000
m = np.zeros(n)
m[:] = 1
m[n-1] = .1 ##benda ke3


t_ev = np.arange(0,t_end,dt)
conv = 1./2*np.pi # year
t_ev_yr = t_ev*conv
# print(t)

## Filename for saved video name
filename = ""

r = np.zeros([n,3])
v = np.zeros([n,3])
a = np.zeros([n,3])
jk = np.zeros([n,3])

def circle():
	filename = "circle-{}-body".format(n)
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
	return filename

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
	return filename


def inf():
	filename = "inf-shape-3-body"
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
	return filename

def plummer():
	filename = "Plumer-{}-body".format(n)
	for i in range(n):
		rad = 1. / np.sqrt( np.random.rand() ** (-2.0/3.0) - 1.0)
		theta = np.arccos(np.random.uniform(-1,1))
		phi = np.random.uniform(0,2*np.pi)
		r[i,0] = rad * np.sin( theta ) * np.cos( phi )
		r[i,1] = rad * np.sin( theta ) * np.sin( phi )
		r[i,2] = rad * np.cos( theta )

	v_abs = np.random.rand() * np.sqrt(2.0) * ( 1.0 + rad*rad)**(-0.25)
	for i in range(n):
		phi = i * 2 * np.pi / n
		v[i,0] = v_abs * np.sin( theta ) * np.cos( phi )
		v[i,1] = v_abs * np.sin( theta ) * np.sin( phi )
		v[i,2] = v_abs * np.cos( theta )
	return filename

def benda_ke3():
	filename = "benda-ke-3"
	for i in range(n):
		phi = i * 2 * np.pi / n
		r[i,0] = np.cos(phi)
		r[i,1] = np.sin(phi)
		r[i,2] = 0
	v_abs = .8/np.sqrt(np.sqrt(n))
	for i in range(n):
		phi = i * 2 * np.pi / n
		v[i,0] = - v_abs * np.sin(phi)
		v[i,1] = v_abs * np.cos(phi)
		v[i,2] = 0
	# Custom this to modify the 3rd body 
	r[n-1,0]=10
	r[n-1,1]=10
	r[n-1,2]=0
	v[n-1,0]=0
	v[n-1,1]=-0.2
	v[n-1,2]=.0
	return filename

filename = benda_ke3()
# circle()

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

# v[0,0] += 0.0001
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
t_out = 0.

old_r=np.zeros([n,3])
old_v=np.zeros([n,3])
old_a=np.zeros([n,3])
old_j=np.zeros([n,3])

rx = []
ry = []
rz = []
tt = []

for t in t_ev:
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
	# fig = plt.figure(0)
	# l, = plt.plot(r[:,0],r[:,1],"b.")

	if t % 1 == 0:		# Saving Frame each 10 step
		rx.append(np.array(r[:,0]))
		ry.append(np.array(r[:,1]))
		rz.append(np.array(r[:,2]))
		tt.append(t_ev_yr[:])
	t += 1

rx = np.array(rx)
ry = np.array(ry)
rz = np.array(rz)
tt = np.array(tt)
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



def update_line(num, rx, ry, line1, line2, tt):
    xx = rx[num, :]
    yy = ry[num, :]
    tev = tt[num, :]
    trailx.append(xx)
    traily.append(yy)
    line1.set_data(xx,yy)
    for i in range(num):
	    time_ev.set_text("t = {:.4f} yr".format(tev[i]))
    line2.set_data(trailx, traily)
    return line1, line2

def update_line_notrail(num, rx, ry, line, tt):
    xx = rx[num, :]
    yy = ry[num, :]
    tev = tt[num, :]
    line.set_data(xx,yy)
    for i in range(num):
	    time_ev.set_text("t = {:.4f} yr".format(tev[i]))
    return line
    
trailx = []
traily = []

fig = plt.figure(0)
ax = fig.add_subplot(111)
l2, = ax.plot([],[],"g.", markersize='1')
l1, = ax.plot([],[],"r*", markersize='5')
time_ev = ax.text(0.8, 0.95, '', transform=ax.transAxes)

plt.xlim(np.min(rx[:,1]),np.max(rx[:,1]))
plt.ylim(np.min(ry[:,1]),1.5*np.max(ry[:,1]))
## batas-batas untuk kluster
# candidate = [np.max(rx[:,1]),abs(np.min(rx[:,1])),np.max(ry[:,1]),abs(np.min(ry[:,1]))]
# lim = np.max(candidate)
# plt.xlim(-1*lim,lim)
# plt.ylim(-1*lim,lim)

plt.xlabel("X (AU)")
plt.ylabel("Y (AU)")
plt.title(filename)

line_ani = animation.FuncAnimation(fig, update_line, len(rx), fargs=(rx,ry,l1,l2,tt), interval=10)
# line_ani = animation.FuncAnimation(fig, update_line_notrail, len(rx), fargs=(rx,ry,l1,tt), interval=10)

e_out = ekin + epot
print ("Final total energy E_out = {} \n".format(e_out))
print ("absolute energy error: E_out - E_in = {} \n".format(e_out - e_in))
print ("relative energy error: (E_out - E_in) / E_in = {} ".format((e_out - e_in) / e_in ))
line_ani.save("{}.mp4".format(filename), fps=60)
# line_ani.save("{}-notrail.mp4".format(filename), fps=60)
print("Finish in %s seconds" % (time.time() - start_time))
# print(t_ev_yr)
plt.show()
