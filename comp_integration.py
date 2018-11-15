"""
I took the prefactor 1/(1+z) and multiplied it to
the result obtained from integrating the integrand.
"""

#Import libraries
from math import sqrt
import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from gaussxw import gaussxw

#Creating lists for graphing
trap_list = []
simp_list = []
gauss_list = []
trapz_list = []
axis = []

#Creating list for numpy.trapz y values
y_list = []

#Defining iteration from 0 to 10
z = -1
while(z < 10):
    z += 1

    #Appending z for graphing on the x-axis
    axis.append(z)
    #Defining function
    def f(w):
        return ((3e8)/(70000*(sqrt(0.3*((1+w)**3) + 0.7))))
    #Appending list of f(z) for numpy.trapz
    y_list.append(f(z))
    #Trapezoid method
    def trapezoid(N,a,b):
        h = (b-a)/N
        s = 0.5*(f(a) + f(b))
        for k in range(1,N):
            s += f(a + k*h)
        return (1/(1+z))*h*s
    result=trapezoid(1000,0,z)
    trap_list.append(result)
    print("\nIntegral with trapezoid rule for z = " ,z, "is\n",result,)

    #Simpson's method
    def simpson(N,a,b):

        h = (b-a)/N
        s = (f(a) + f(b))

        for k in range(1,N,2):
            s += 4 * f(a + k*h)
        for k in range(2,N,2):
            s += 2 * f(a + k*h)

        return (1/(1+z))*s*h/3

    result=simpson(1000,0,z)
    simp_list.append(result)
    print("\nIntegral with simpson's rule for z =",z,":\n",result,)

    #Gaussian quadrature method
    def gaussq(N,a,b):

        #Calculation of points and weights
        x,w = gaussxw(N)
        xp = 0.5 * (b-a)*x + 0.5*(b+a)
        wp = 0.5 *(b-a)*w

        #Integration
        s = 0.0
        for k in range(N):
            s += wp[k] * f(xp[k])

        return (1/(1+z))*s


    result=gaussq(1000,0,z)
    gauss_list.append(result)
    print("\nIntegral with gaussian quadrature for N=1000 and z =",z," :\n",result,)

    #np.trapz method
    def nptrapz(y, a, dz):
        return (1/(1+z))*np.trapz(y, x=a, dx=dz)
    result=nptrapz(y_list, axis, 1)
    print("\nIntegral with np.trapz method for z =",z," :\n",result,)
    trapz_list.append(result)

#Graphing functions
def graph_trap(func1, x1_range):
    x = np.array(x1_range)
    y = np.array(func1)
    fig, ax = plt.subplots()
    pl.ylim([800,1800])
    pl.xlim([0,10])
    plt.margins(0)
    ax.plot(x, y)
    ax.set(xlabel='$z$', ylabel='$D_A(z)$ [Mpc]',title='Trapezoid Rule')
    plt.savefig('graph_trap.pdf')
    plt.show()


graph_trap(trap_list, axis)

def graph_simp(func2, x2_range):
    x = np.array(x2_range)
    y = np.array(func2)
    fig, ax = plt.subplots()
    pl.ylim([800,1800])
    plt.margins(0)
    ax.plot(x, y)
    ax.set(xlabel='$z$', ylabel='$D_A(z)$ [Mpc]',title="Simpson's Rule")
    plt.savefig('graph_simp.pdf')
    plt.show()

graph_simp(simp_list, axis)

def graph_gauss(func3, x3_range):
    x = np.array(x3_range)
    y = np.array(func3)
    fig, ax = plt.subplots()
    pl.ylim([800,1800])
    plt.margins(0)
    ax.plot(x, y)
    ax.set(xlabel='$z$', ylabel='$D_A(z)$ [Mpc]',title='Gaussian Quadrature Rule')
    plt.savefig('graph_gauss.pdf')
    plt.show()

graph_gauss(gauss_list, axis)


def graph_trapz(func4, x4_range):
    x = np.array(x4_range)
    y = np.array(func4)
    fig, ax = plt.subplots()
    pl.ylim([800,1800])
    plt.margins(0)
    ax.plot(x, y)
    ax.set(xlabel='$z$', ylabel='$D_A(z)$ [Mpc]',title='numpy.trapz Rule')
    plt.savefig('graph_trapz.pdf')
    plt.show()

graph_trapz(trapz_list, axis)
