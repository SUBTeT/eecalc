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


# main
if __name__ == '__main__':
    A = np.array([[1.0, 2.0, 0.0], [1.0, -1.0, 1.0], [2.0, 0.0, 1.0]])
    b = np.array([1, -1, -2])

    y = calc_y(lu(A), b)
    print("A:\n", A)
    print("b:\n", b)
    print("y:\n", y)
