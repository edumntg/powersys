class RK4(object):
    #4th-order Runge-Kutta method

    @staticmethod
    def step(f, x, y, h, *argv):
        """Solves a system of ordinary differential equations using the 4th-order Runge-Kutta method.

        Args:
            f: A function that calculates the derivatives of the state variables.
            x: The independent variable.
            y: A NumPy array containing the current values of the state variables.
            h: The time step (seconds).
            *argv: A list containing all other arguments (optinal).

        Returns:
            A NumPy array containing the new values of the state variables.
        """
        # Compute Rk4 terms
        k1 = f(x, y, *argv)
        k2 = f(x + 0.5*h, y + 0.5*h*k1, *argv)
        k3 = f(x + 0.5*h, y + 0.5*h*k2, *argv)
        k4 = f(x + h, y + h*k3, *argv)

        return y + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)