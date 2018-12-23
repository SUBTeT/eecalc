import numpy as np
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


def solve_by_jacobi(A: np.ndarray, b: np.ndarray, x: np.ndarray) -> np.ndarray:
    while True:
        x_next = calc_next(A, b, x)
        dif = np.abs(np.linalg.norm(x_next) - np.linalg.norm(x))
        if dif <= error:
            break
        if dif > limit:
            print("Jacobi method cannot get the solution...")
            exit(1)
        x = x_next
    return x_next


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
    # size of matrix
    size = 5
    a = -0.25
    A, b = make(size, a)
    x = np.zeros(size)
    print(solve_by_jacobi(A, b, x))
