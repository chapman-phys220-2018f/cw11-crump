#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
# Name: Conner Carnahan
# Student ID: 1614309
# Email: carna104@mail.chapman.edu
# Course: PHYS220/MATH220/CPSC220 Fall 2018
# Assignment: HW 10
###

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


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
    
    uk = npt
    upk = npt
        
    count = 0
    uk[0] = f0
    upk[0] = f1
    for ti in np.nditer(npt):
        uk[count] = uk[count-1]+dt*upk[count-1]
        upk[count] = upk[count-1]-dt*uk[count-1]
        count += 1
        
    plotitokay(npt,uk,upk)

    
def plotitokay(t,uk,upk):
    fig = plt.figure(figsize = (12,8))
    a = plt.axes()
    
    a.plot(t,uk, 'b', label = "$u(t)$")
    a.plot(t,upk, 'r', label = "$u'(t)$")
    
    a.legend()
    plt.show()