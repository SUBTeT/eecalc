import sympy as sp
import copy as cp

# equation
def p_0(x):
    return 1
def p_1(x):
    return x

y = sp.Symbol('y')

# Array of Legendre polynomial : p_n
p_n = [p_0(y), p_1(y)]
for i in range(1, 10):
    p_i = ((2*i + 1)*y*p_n[-1] - i*p_n[-2])/(i+1)
    p_n.append(sp.simplify(p_i))

# Dichotomy
def dichotomy(func_f, x_min, x_max, i):
    num_calc = 0  # the number of calculations

    while(True):
        # calculate new x_mid
        x_mid = (x_max+x_min)/2.0

        # refresh the search interval
        if (0.0 < func_f.subs([(y, x_mid)])*func_f.subs([(y, x_max)])):  # if sign of middle == right
            x_max = x_mid  #refresh the right
        else:  # if sign of middle == left
            x_min = x_mid  # refresh the left

        num_calc += 1

        # end if the number of calculations is greater than i
        if(i<=num_calc):
            break

    return x_mid # return the answer

# solve the Legendre polynomial: P[number] with dichotomy
def calculate(number, a):
    b = []
    if number%2==1: # n is odd
        b.append(0)
    for i in range(0, number//2-1):
        sol = dichotomy(p_n[number], a[i], a[i+1], 30)
        b.append(sol)
    sol = dichotomy(p_n[number], a[number//2-1], 1, 30)
    b.append(sol)
    return b

def show(b, i): # manipulation of arrays
    c = cp.copy(b)

    for j in range(len(c)):
        if c[j] != 0:
            c.append(-c[j])
    c.sort()
    if i >= 4:
        print("P[{}] : {}".format(i, c))

    return b

# main
if __name__ == '__main__':
    a = [0] # initial setting

    print("Zero points of Legendre polynomial P[n]")
    for i in range(2, 11):
        b = calculate(i, a)
        a = show(b, i)


