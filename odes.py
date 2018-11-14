#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
# Name: Conner Carnahan
# Student ID: 1614309
# Email: carna104@mail.chapman.edu
#FULL NAME :NATANAEL ALPAY
#ID        :002285534
#email:alpay100@mail.chapman.edu

# Course: PHYS220/MATH220/CPSC220 Fall 2018
# Assignment: HW 10
###

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
import array as ar


def noteulermethod(f,g,t0,tf,f0,f1,dt = 0.001):
    """eulermethod(sympy function f, sympy function g, float t0, float tf, float f0, float f1, float dt = 0.001): 
       returns: numpy array of two arrays
       does an euler approximation for a system of differential equations u = f[u,u',t], u' = g[u,u',t]
       """
    u,up,t = sp.symbols('u up t')
    
    npt = np.linspace(float(t0),float(tf),int((tf-t0)/dt))
    
    uk = npt
    upk = npt
    
    npf = sp.lambdify(t,u,up,f)
    npg = sp.lambdify(t,u,up,g)
    
    count = 0
    uk[0] = f0
    upk[0] = f1
    for ti in np.nditer(npt):
        uk[count] = npf(uk[count-1],upk[count-1],float(ti))
        upk[count] = npg(uk[count-1],upk[count-1],float(ti))
        count += 1
        
    plotitokay(npt,uk,upk)

def eulermethod(t0,tf,f0,f1,dt = 0.001):
    """eulermethod(sympy function f, sympy function g, float t0, float tf, float f0, float f1, float dt = 0.001): 
       returns: numpy array of two arrays
       does an euler approximation for a system of differential equations u = f[u,u',t], u' = g[u,u',t]
       """

    npt = np.linspace(float(t0),float(tf),int((tf-t0)/dt))
    
    uk = np.zeros(int((tf-t0)/dt))
    upk = np.zeros(int((tf-t0)/dt))
        
    count = 1
    uk[0] = float(f0)
    upk[0] = float(f1)
    
    while (count < npt.size):
        uk[count] = uk[count-1] + dt*upk[count-1]
        upk[count] = upk[count-1] - dt*uk[count-1]
        count += 1
    
    plotitokay(npt,uk,upk)


def heun(t0,tf,f0,f1,dt = 0.001):
    
    npt = np.linspace(float(t0),float(tf),int((tf-t0)/dt))
    
    uk = np.zeros(int((tf-t0)/dt))
    upk = np.zeros(int((tf-t0)/dt))
    
    ubar = 0.0
    upbar = 0.0
    count = 1
    
    uk[0] = f0
    upk[0] = f1
    
    while (count < npt.size):
        ubar = uk[count-1] + dt*upk[count-1]
        upbar = upk[count-1] - dt*uk[count-1]
        uk[count] = uk[count-1] + (dt/2)*(upk[count-1]+ubar)
        upk[count] = uk[count-1] - (dt/2)*(upbar+uk[count-1])
        count += 1
        
    plotitokay(npt,uk,upk)
    
    
    
#def RungeKutta_2(0,tf,f0,f1,dt = 0.001):
def RungeKutta_2( f, t0, t ):
    n = len( t )
    rk = np.ar([ t0 ]*n)
    while count<n-1:
        h = t[count+1] - t[count]
        k1 = h * f( x[count], t[count] ) / 2.0
        rk[count+1] = x[count] + h * f( x[count] + k1, t[count] + h / 2.0 )
        count=count+1
    """
        f     - function 
        x0    - the initial condition
        t     - list or NumPy array of t values to compute solution at.

    """
def RungeKutta_4( f, t0, t ):
    count=0
    n = len(t)
    rk4 = np.ar([t0]* n)
    while count<n-1:
        h = t[count+1] - t[count]
        k1 = h * f( rk4[count], t[i] )
        k2 = h * f( rk4[count] + 0.5 * k1, t[count] + 0.5 * h )
        k3 = h * f( rk4[count] + 0.5 * k2, t[count] + 0.5 * h )
        k4 = h * f( rk4[count] + k3, t[count+1] )
        rk4[count+1] = rk4[count] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0
        count=count+1
        

    
def plotitokay(t,uk,upk):
    fig = plt.figure(figsize = (16,9))
    a = plt.axes()
    
    a.plot(t,uk, 'b', label = "$u(t)$")
    a.plot(t,upk, 'r', label = "$u'(t)$")
    
    a.legend()
    plt.show()