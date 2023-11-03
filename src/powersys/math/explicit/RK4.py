class RK4(object):
    #4th-order Runge-Kutta method

    @staticmethod
    def step(f, x, y, h):
        # Compute Rk4 terms
        k1 = f(x, y)
        k2 = f(x + 0.5*h, y + 0.5*h*k1)
        k3 = f(x + 0.5*h, y + 0.5*h*k2)
        k4 = f(x + h, y + h*k3)

        return y + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)