import numpy as np


####################################################################################################
# Exercise 1: Interpolation

def lagrange_interpolation(x: np.ndarray, y: np.ndarray) -> (np.poly1d, list):
    """
    Generate Lagrange interpolation polynomial.

    Arguments:
    x: x-values of interpolation points
    y: y-values of interpolation points

    Return:
    polynomial: polynomial as np.poly1d object
    base_functions: list of base polynomials
    """

    assert (x.size == y.size)

    polynomial = np.poly1d(0)
    base_functions = []
    # TODO: Generate Lagrange base polynomials and interpolation polynomial

    # Generate Lagrange base polynomials and interpolation polynomial
    n = x.size  # number of points
    for i in range(n):
        # create a base polynomial of degree n-1
        base = np.poly1d(1)
        for j in range(n):
            if i != j:
                # multiply by a factor of (x - x_j) / (x_i - x_j)
                base *= np.poly1d([1, -x[j]]) / (x[i] - x[j])
        # add the scaled base polynomial to the interpolation polynomial
        polynomial += y[i] * base
        # append the base polynomial to the list
        base_functions.append(base)

    return polynomial, base_functions




def hermite_cubic_interpolation(x: np.ndarray, y: np.ndarray, yp: np.ndarray) -> list:
    """
    Compute Hermite cubic interpolation spline.

    Arguments:
    x (np.ndarray): x-values of interpolation points
    y (np.ndarray): y-values of interpolation points
    yp (np.ndarray): derivative values of interpolation points

    Returns:
    List[np.poly1d]: List of np.poly1d objects, each interpolating the function between two adjacent points
    """
    assert (len(x) == len(y) == len(yp)), "x, y, and yp arrays must have the same length."

    n = len(x) - 1
    spline = []

    for i in range(n):
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]

        # Matrix form to solve for coefficients
        M = np.array([
            [1, x[i], x[i]**2, x[i]**3],
            [1, x[i+1], x[i+1]**2, x[i+1]**3],
            [0, 1, 2*x[i], 3*x[i]**2],
            [0, 1, 2*x[i+1], 3*x[i+1]**2]
        ])

        b = np.array([y[i], y[i+1], yp[i], yp[i+1]])

        # Solve for coefficients [d, c, b, a]
        coeffs = np.linalg.solve(M, b)

        # Create the cubic polynomial for this segment
        poly = np.poly1d(coeffs[::-1])
        spline.append(poly)

    return spline




####################################################################################################
# Exercise 2: Animation

def natural_cubic_interpolation(x: np.ndarray, y: np.ndarray) -> list:
    """
    Intepolate the given function using a spline with natural boundary conditions.

    Arguments:
    x: x-values of interpolation points
    y: y-values of interpolation points

    Return:
    spline: list of np.poly1d objects, each interpolating the function between two adjacent points
    """
    assert (x.size == y.size)

    # TODO construct linear system with natural boundary conditions

    # TODO solve linear system for the coefficients of the spline

    n = len(x) - 1
    h = np.diff(x)
    b = np.diff(y) / h

    # Set up the linear system for the second derivatives
    A = np.zeros((n+1, n+1))
    v = np.zeros(n+1)

    # Filling the A matrix
    A[0, 0], A[n, n] = 1, 1
    for i in range(1, n):
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
        v[i] = 6 * (b[i] - b[i-1])

    # Solve for the second derivatives
    m = np.linalg.solve(A, v)

    # Compute spline coefficients
    spline = []
    # TODO extract local interpolation coefficients from solution
    for i in range(n):
        a = (m[i+1] - m[i]) / (6 * h[i])
        b = m[i] / 2
        c = -(h[i] * m[i+1]) / 6 - (h[i] * m[i]) / 3 + (y[i+1] - y[i]) / h[i]
        d = y[i]

        spline.append(np.poly1d([a, b, c, d]))

    return spline



def periodic_cubic_interpolation(x: np.ndarray, y: np.ndarray) -> list:
    """
    Interpolate the given function with a cubic spline and periodic boundary conditions.

    Arguments:
    x: x-values of interpolation points
    y: y-values of interpolation points

    Return:
    spline: list of np.poly1d objects, each interpolating the function between two adjacent points
    """
    assert (x.size == y.size)
    # TODO: construct linear system with periodic boundary conditions

    # TODO solve linear system for the coefficients of the spline

    n = len(x)
    h = np.diff(x)
    b = np.diff(y) / h

    # Adjust for periodicity
    h = np.append(h, h[0])
    b = np.append(b, (y[0] - y[-1]) / h[-1])

    # Set up the linear system for the second derivatives
    A = np.zeros((n, n))
    v = np.zeros(n)

    A[0, 0], A[0, -1], A[0, 1] = 2*(h[-1] + h[0]), h[-1], h[0]
    A[-1, 0], A[-1, -2], A[-1, -1] = h[0], h[-2], 2*(h[-2] + h[-1])
    for i in range(1, n-1):
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
        v[i] = 6 * (b[i] - b[i-1])
    v[0] = 6 * (b[0] - b[-1])
    v[-1] = v[0]

    # Solve for the second derivatives
    m = np.linalg.solve(A, v)

    # Compute spline coefficients
    spline = []
    # TODO extract local interpolation coefficients from solution
    for i in range(n-1):
        a = (m[i+1] - m[i]) / (6 * h[i])
        b = m[i] / 2
        c = -(h[i] * m[i+1]) / 6 - (h[i] * m[i]) / 3 + (y[i+1] - y[i]) / h[i]
        d = y[i]

        spline.append(np.poly1d([a, b, c, d]))

    return spline



if __name__ == '__main__':

    x = np.array( [1.0, 2.0, 3.0, 4.0])
    y = np.array( [3.0, 2.0, 4.0, 1.0])

    splines = natural_cubic_interpolation( x, y)

    # # x-values to be interpolated
    # keytimes = np.linspace(0, 200, 11)
    # # y-values to be interpolated
    # keyframes = [np.array([0., -0.05, -0.2, -0.2, 0.2, -0.2, 0.25, -0.3, 0.3, 0.1, 0.2]),
    #              np.array([0., 0.0, 0.2, -0.1, -0.2, -0.1, 0.1, 0.1, 0.2, -0.3, 0.3])] * 5
    # keyframes.append(keyframes[0])
    # splines = []
    # for i in range(11):  # Iterate over all animated parts
    #     x = keytimes
    #     y = np.array([keyframes[k][i] for k in range(11)])
    #     spline = natural_cubic_interpolation(x, y)
    #     if len(spline) == 0:
    #         animate(keytimes, keyframes, linear_animation(keytimes, keyframes))
    #         self.fail("Natural cubic interpolation not implemented.")
    #     splines.append(spline)

    print("All requested functions for the assignment have to be implemented in this file and uploaded to the "
          "server for the grading.\nTo test your implemented functions you can "
          "implement/run tests in the file tests.py (> python3 -v test.py [Tests.<test_function>]).")
