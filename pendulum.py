# Goal: write the simplest explicit Python program to solve the pendulum probelm
import pylab as pl

def acc(x,g,l):
  """Pendulum'sangular acceleration"""
  return (-g/l)*pl.sin(x)

def update(x,v,a,g,l,dt):
  """Velocity Verlet integrator"""
  x += v*dt + 0.5 *a*dt**2
  v += 0.5* (acc(x, g, l)+acc(x, g, l))*dt
  return x,v,a

def analytic(t, g,l, x0, v0):
  """Small angle analytic solution"""
  omega = pl.sqrt(g/l)
  return x0*pl.cos(omega*t)+v0/omega*pl.sin(omega*t)

g =9.81
l = 1.0
dt = 0.001
totaltime = 20
nsteps = int(totaltime/dt)
t = pl.arange(nsteps)*dt

# initial condition
x,v=x0,v0=0.1,0.0
a = acc(x, g, l)
X ,V,A=[],[],[]
for i in range(nsteps):
  X.append(x), V.append(v), A.append(a)
  x,v,a = update(x, v, a, g, l, dt)
  

pl.plot(t,X,'k-', lw=3, label="numerical")
pl.plot(t,analytic(t,g,l,x0,v0), lw=1, label="small angle")
pl.xlabel("time"), pl.ylabel(r"$\theta$"), pl.legend(loc="right", frameon=True, fontsize="x-small")
pl.show()