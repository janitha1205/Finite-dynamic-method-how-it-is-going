import matplotlib.pyplot as plt
import numpy as np


class element:

    g = np.array([0, 9.81, 0]).T

    k11 = 0.103
    k12 = 0.002
    k13 = 0.004
    k22 = 0.109
    k23 = 0.007
    k33 = 0.103
    b11 = 0.103
    b12 = 0.002
    b13 = 0.004
    b22 = 0.109
    b23 = 0.007
    b33 = 0.103
    k = np.array([[k11, k12, k13], [k12, k22, k23], [k13, k23, k33]])
    b = np.array([[b11, b12, b13], [b12, b22, b23], [b13, b23, b33]])
    dt = 0.01

    def __init__(self, Fx, Fy, Fz, m):
        self.F = np.array([Fx, Fy, Fz]).T
        self.xdot = np.array([0, 0, 0]).T
        self.x = np.array([0, 0, 0]).T

        self.x2dot = np.array([0, 0, 0]).T
        self.x2 = np.array([0, 0, 0]).T
        self.m = m

    def cal_df(self):
        xdotdot = (
            -np.dot(self.k, self.x)
            - np.dot(self.b, self.xdot)
            + self.F
            - self.g * self.m
        ) * (2 / self.m)
        print(xdotdot.shape)
        self.xdot = xdotdot * self.dt + self.xdot
        self.x = 0.5 * xdotdot * self.dt * self.dt + self.xdot * self.dt + self.x
        dF = (
            self.F
            - (self.m / 2) * xdotdot
            - np.dot(2 * self.k, self.x)
            - np.dot(self.b / 2, self.xdot)
        )
        return dF

    def cal_x_2(self, dF):
        x2dotdot = (
            dF - np.dot(2 * self.k, self.x2) - np.dot(self.b / 2, self.x2dot)
        ) / (self.m / 2)
        self.x2dot = x2dotdot * self.dt + self.x2dot
        self.x2 = 0.5 * x2dotdot * self.dt * self.dt + self.x2dot * self.dt + self.x2

    def run_simulation(self, T):
        t = 0
        tt = []
        x11 = []
        x21 = []
        while t < T:
            df = self.cal_df()
            self.F = np.array([0, 0, 0]).T
            self.cal_x_2(df)
            tt.append(t)
            x11.append(self.x)
            x21.append(self.x2)
            t = t + self.dt
        return tt, x11, x21


el1 = element(1.56, 1.3, 1.7, 0.02)
t, x1, x2 = el1.run_simulation(10)
plt.plot(t, x1, t, x2)


plt.show()
