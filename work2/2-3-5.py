import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numba import jit


# LU decomposition
@jit
def lu(A: np.ndarray) -> np.ndarray:
    n = len(A[0])
    LU = np.zeros(shape=(n, n))
    for i in range(0, n):
        LU[0][i] = A[0][i]
    for i in range(1, n):
        LU[i][0] = A[i][0] / LU[0][0]

    for j in range(1, n):
        for i in range(1, j + 1):
            sum = 0
            for k in range(0, i):
                sum += LU[i][k] * LU[k][j]
            LU[i][j] = A[i][j] - sum
        for i in range(j + 1, n):
            sum = 0
            for k in range(0, j):
                sum += LU[i][k] * LU[k][j]
            LU[i][j] = (A[i][j] - sum) / LU[j][j]
    return LU


@jit
def calc_y(LU: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = len(b)
    y = np.zeros(n)
    y[0] = b[0]
    for i in range(1, n):
        sum = 0
        for k in range(0, i):
            sum += LU[i][k] * y[k]
        y[i] = b[i] - sum
    return y


@jit
def calc_x(LU: np.ndarray, y: np.ndarray) -> np.ndarray:
    n = len(y)
    x = np.zeros(n)
    x[n - 1] = y[n - 1] / LU[n - 1][n - 1]
    for i in reversed(range(0, n - 1)):
        sum = 0
        for k in reversed(range(i + 1, n)):
            sum += LU[i][k] * x[k]
        x[i] = (y[i] - sum) / LU[i][i]
    return x


def solve_by_lu(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    LU = lu(A)
    y = calc_y(LU, b)
    x = calc_x(LU, y)
    return x


def make_plane(g: int, w: int, d: int, l: int, h: float,
               x: sp.symbol.Symbol) -> dict:  # set the "coordinates : value -- dictionary"
    plane = {}
    for i in range(0, int(l // (2 * h) + 1)):  # A
        plane[(i, int(g // h))] = 100

    for i in range(int((w / 2) // h), int((l / 2) // h + 1)):  # B
        plane[(i, 0)] = 0
    for i in range(0, int((w / 2) // h) + 1):
        plane[(i, int(-d // h))] = 0
    for j in range(0, int(d // h + 1)):
        plane[((w / 2) // h, -j)] = 0

    for j in range(1, int(g // h)):  # C
        plane[(int((l / 2) // h), j)] = j * h * 100 / g

    k = 1
    for i in range(0, int((w / 2) // h)):  # D
        for j in range(int(-d // h + 1), int(g // h)):
            plane[(i, j)] = x ** k
            k += 1

    for j in range(1, int(g // h)):  # E
        for i in range(int((w / 2) // h), int((l / 2) // h)):
            plane[(i, j)] = x ** k
            k += 1

    for j in range(int(-d // h), int(g // h + 1)):  # F
        plane[(-1, j)] = plane[(1, j)]

    return plane


def make_f_list(g: int, w: int, d: int, l: int, h: float, plane: dict) -> list:  # make the equation list
    f_list = []  # the functions are 4 * center - (surroundings)

    for i in range(0, int((w / 2) // h)):  # D
        for j in range(int(-d // h + 1), int(g // h)):
            f_list.append(
                4 * plane[(i, j)] - plane[(i + 1, j)] - plane[(i - 1, j)] - plane[(i, j + 1)] - plane[(i, j - 1)])

    for j in range(1, int(g // h)):  # E
        for i in range(int((w / 2) // h), int((l / 2) // h)):
            f_list.append(
                4 * plane[(i, j)] - plane[(i + 1, j)] - plane[(i - 1, j)] - plane[(i, j + 1)] - plane[(i, j - 1)])

    return f_list


def make_matrix(g: int, w: int, d: int, l: int, h: float) -> (
        np.ndarray, np.ndarray, dict):  # make the coefficient matrix
    x = sp.Symbol('x')
    plane = make_plane(g, w, d, l, h, x)  # plane is the "coordinates : value -- dictionary"
    f_list = make_f_list(g, w, d, l, h, plane)
    n = len(f_list)  # number of variables
    A = []
    b = []

    for f in f_list:
        line = sp.poly(f).all_coeffs()  # fetch the coefficients
        b.append(-line.pop())  # pop the last elements into b
        line.reverse()  # adjust order of the list

        if n - len(line) != 0:  # adjust the line size to n
            for i in range(n - len(line)):
                line.append(0)

        A.append(line)  # append line to A

    return np.array(A), np.array(b), plane


def plot_plane(g: int, w: int, d: int, l: int, h: float, plane: dict, sol: np.ndarray):
    del_lst = []
    for tuple, value in plane.items():  # delete (-1, j)
        if tuple[0] == -1:
            del_lst.append(tuple)

    for tuple in del_lst:
        del plane[tuple]

    k = 0
    for i in range(0, int((w / 2) // h)):  # D
        for j in range(int(-d // h + 1), int(g // h)):
            plane[(i, j)] = sol[k]
            k += 1

    for j in range(1, int(g // h)):  # E
        for i in range(int((w / 2) // h), int((l / 2) // h)):
            plane[(i, j)] = sol[k]
            k += 1
    x = []
    y = []
    z = []
    for tuple, value in plane.items():
        x.append(tuple[0])
        y.append(tuple[1])
        z.append(value)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_trisurf(x, y, z, cmap='plasma')
    plt.show()


if __name__ == '__main__':
    g, w, d, l, h = 10, 10, 10, 30, 5
    A, b, plane = make_matrix(g, w, d, l, h)

    x = solve_by_lu(A, b)
    print(x)
    plot_plane(g, w, d, l, h, plane, x)
