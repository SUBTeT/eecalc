import numpy as np
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


# make the matrices A, b
def make(n: int, a: float) -> (np.ndarray, np.ndarray):
    A = np.zeros(shape=(n, n))
    b = np.zeros(n)

    b[0] = 1.
    b[n - 1] = -1.

    A[0][n - 1] = a
    A[n - 1][0] = a
    for i in range(n):
        A[i][i] = 1.
    for i in range(n - 1):
        A[i][i + 1] = a
    for i in range(1, n):
        A[i][i - 1] = a
    return A, b


# main
if __name__ == '__main__':
    # size of matrix
    size = 5
    a = -0.75
    A, b = make(size, a)
    x = solve_by_lu(A, b)
    print(x)
