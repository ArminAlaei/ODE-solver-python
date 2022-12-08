

from ODESolver import *
import matplotlib.pyplot as plt
class Heuns(ODESolver):
    def advance(self):
        u, f, k, t = self.u, self.f, self.k, self.t

        dt =  t[k+1] - t[k]
        k1 = f(u[k],t[k])
        k2 = f(u[k] + dt*k1, t[k] + dt)
        unew = u[k] + dt*((k1/2) + (k2/2))
        return unew

class Midpoint(ODESolver):
    def advance(self):
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k+1] - t[k]
        dt2 = dt/2.0
        k1 = f(u[k], t)
        k2 = f(u[k] + dt2*k1, t[k] + dt2)
        unew = u[k] + dt*k2
        return unew


#defining the function
def f(u,t):
    return t * np.cos(t) - np.sin(t)



t1 = np.linspace(0 , np.pi * 8, 20)
t2 = np.linspace(0 , np.pi * 8, 25)
t3 = np.linspace(0 , np.pi * 8, 50)
t4 = np.linspace(0 , np.pi * 8, 150)

#exact solution
ex =  t1 * np.sin(t1) -  2*np.cos(t1)

#defining different methods
H = Heuns(f)
H.set_initial_condition(2)

Mid = Midpoint(f)
Mid.set_initial_condition(2)

rk4 = RungeKutta4(f)
rk4.set_initial_condition(2)

#Solving using different methods
uH1, tH1 = H.solve(t1)
uH2, tH2 = H.solve(t2)
uH3, tH3 = H.solve(t3)
uH4, tH4 = H.solve(t4)

uM1, tM1 = Mid.solve(t1)
uM2, tM2 = Mid.solve(t2)
uM3, tM3 = Mid.solve(t3)
uM4, tM4 = Mid.solve(t4)

uR1, tR1 = rk4.solve(t1)
uR2, tR2 = rk4.solve(t2)
uR3, tR3 = rk4.solve(t3)
uR4, tR4 = rk4.solve(t4)

#Plotting for different methods
plt.plot(tH1,uH1, label= "n = 20")
plt.plot(tH2,uH2, label= "n = 25")
plt.plot(tH3,uH3, label= "n = 50")
plt.plot(tH4,uH4, label= "n = 150")
plt.plot(t1,ex, label = "exact")
plt.grid()
plt.title("Heuns methods")
plt.legend()
plt.show()

plt.plot(tM1,uM1, label= "n = 20")
plt.plot(tM2,uM2, label= "n = 25")
plt.plot(tM3,uM3, label= "n = 50")
plt.plot(tM4,uM4, label= "n = 150")
plt.plot(t1,ex, label = "exact")
plt.grid()
plt.title("Midpoint method")
plt.legend()
plt.show()

plt.plot(tR1,uR1, label= "n = 20")
plt.plot(tR2,uR2, label= "n = 25")
plt.plot(tR3,uR3, label= "n = 50")
plt.plot(tR4,uR4, label= "n = 150")
plt.plot(t1,ex, label = "exact")
plt.grid()
plt.title("Runge kutta method")
plt.legend()
plt.show()

"""

Run example:

user$ python3 compare_methods.py

plots attached

"""
