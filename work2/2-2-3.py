import numpy as np
import time
import matplotlib.pyplot as plt
from numba import jit

error = 1e-7
limit = 1e+7


@jit
def calc_next(A: np.ndarray, b: np.ndarray, x: np.ndarray) -> np.ndarray:
    n = len(b)
    x_next = np.zeros(n)
    for i in range(0, n):
        sum = 0
        for j in range(0, i):
            sum += A[i][j] * x[j]
        for j in range(i + 1, n):
            sum += A[i][j] * x[j]
        x_next[i] = (b[i] - sum) / A[i][i]
    return x_next


def solve_by_jacobi(A: np.ndarray, b: np.ndarray, x: np.ndarray) -> (np.ndarray, float):
    start = time.perf_counter()
    while True:
        x_next = calc_next(A, b, x)
        dif = np.abs(np.linalg.norm(x_next) - np.linalg.norm(x))
        if dif <= error:
            break
        if dif > limit:
            print("Jacobi method cannot get the solution...")
            exit(1)
        x = x_next
    end = time.perf_counter() - start
    return x_next, end


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


def solve_by_lu(A: np.ndarray, b: np.ndarray) -> (np.ndarray, float):
    start = time.perf_counter()
    LU = lu(A)
    y = calc_y(LU, b)
    x = calc_x(LU, y)
    end = time.perf_counter() - start
    return x, end


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


if __name__ == '__main__':
    # input larger than 10
    max_size = 300
    a = -0.25

    time_j = []
    time_lu = []
    x_axis = []
    for i in range(10, max_size + 1):
        x_axis.append(i)

    for N in range(3, max_size + 1):
        A, b = make(N, a)
        x = np.zeros(N)
        x_j, end_j = solve_by_jacobi(A, b, x)
        if (N >= 10):
            time_j.append(end_j)
        x_lu, end_lu = solve_by_lu(A, b)
        if (N >= 10):
            time_lu.append(end_lu)

    plt.scatter(x_axis, time_j, label="jacobi", s=5)  # plot
    plt.scatter(x_axis, time_lu, label="LU", s=5)
    plt.xlabel("n")
    plt.ylabel("time[s]")
    plt.legend()
    plt.show()
