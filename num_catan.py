import numpy as np
import matplotlib.pyplot as plt
import math

c_1 = []
c_2 = []
x = np.arange(1,10)

def catan_formula(n):
    return math.factorial(int(2*n))/((math.factorial(n+1))*(math.factorial(n)))
y = [catan_formula(int(n)) for n in x]

def catan_truco(n):
    s = 1
    for i in range(1,n+1):
        s *= (4*n +2)/(n+2)
    return s
w = [catan_truco(int(n)) for n in x]

print(w)
plt.plot(x, w)
#plt.show()





