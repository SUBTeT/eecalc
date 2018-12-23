import numpy as np
import matplotlib.pyplot as plt

# equation
def func_f(x):
    return (5*x**3-3*x)/2

# dichotomy
def dichotomy(func_f, x_min, x_max, i):
    num_calc = 0  # the number of calculations
    x_tru = np.sqrt(3/5) # true solution
    y = []

    while(True):
        # calculate new x_mid
        x_mid = (x_max+x_min)/2.0
        dif = abs(x_tru-x_mid)

        # refresh the search interval
        if (0.0 < func_f(x_mid)*func_f(x_max)):  # if sign of middle == right
            x_max = x_mid  # refresh the right
        else:  # if sign of middle == left
            x_min = x_mid  # refresh the left

        num_calc += 1

        # end if the number of calculations is greater than i
        if(i<=num_calc):
            break

        y.append(np.log2(dif))

    print("x = {}".format(x_mid))
    plt.xlabel("number of trials")
    plt.ylabel("log_2(error)")
    plt.plot(y) #plot
    plt.show()

    return x_mid # return the answer

# main
if __name__ == '__main__':
    # calculate solution (equation, the left of the interval, the right, the number of calculations）
    solution = dichotomy(func_f, 0.1, 0.9, 40)
