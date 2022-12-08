import numpy as np
import matplotlib.pyplot as plt

class ODESolver:
    def __init__(self, f):
        self.f = f
    def advance(self):
        """Advance solution one time step."""
        raise NotImplementedError # implement in subclass

    def set_initial_condition(self, U0):
        self.U0 = float(U0)

    def solve(self, time_points):
        self.t = np.asarray(time_points)
        N = len(self.t)
        self.u = np.zeros(N)
        # Assume that self.t[0] corresponds to self.U0
        self.u[0] = self.U0
        # Time loop
        for n in range(N-1):
            self.n = n
            self.u[n+1] = self.advance()
        return self.u, self.t

class ForwardEuler(ODESolver):
    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        unew = u[n] + dt*f(u[n], t[n])
        return unew

#defining the function
def f(u, t):
    return u/5



t_vals = [1,2,3,4,5]
fe = ForwardEuler(f)

fe.set_initial_condition(U0=0.1)
for ts in t_vals:
    time_points = np.linspace(0, 20,ts )
    u1, t1 = fe.solve(time_points)
    ex = 0.1*np.exp(t1/5)
    plt.plot(t1,u1, label = "T = "+str(ts))
plt.plot(t1,ex, label = "exact")
plt.legend()
plt.grid()
plt.show()

"""
Run example:

user$ python3 simple_ODE_class_ODESolver.py

plots attached
"""
