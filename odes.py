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

def eulermethod(t0,tf,f0,f1,N):

    npt = np.linspace(float(t0),float(tf),N)
    
    uk = np.zeros(N)
    upk = np.zeros(N)
    dt = float((tf-t0)/N)
    count = 1
    uk[0] = float(f0)
    upk[0] = float(f1)
    
    while (count < npt.size):
        uk[count] = uk[count-1] + dt*upk[count-1]
        upk[count] = upk[count-1] - dt*uk[count-1]
        count += 1
    
    plotitokay(npt,uk,upk, "Euler's Method for N = " + str(N))


def heun(t0,tf,f0,f1,N):
    
    npt = np.linspace(float(t0),float(tf),N)
    
    dt = float((tf-t0)/N)
    uk = np.zeros(N)
    upk = np.zeros(N)
    
    ubar = 0.0
    upbar = 0.0
    count = 1
    
    uk[0] = f0
    upk[0] = f1
    
    while (count < npt.size):
        ubar = uk[count-1] + dt*upk[count-1]
        upbar = upk[count-1] - dt*uk[count-1]
        uk[count] = uk[count-1] + (dt/2)*(upk[count-1] + upbar)
        upk[count] = upk[count-1] - (dt/2)*(ubar + uk[count-1])
        count += 1
        
    plotitokay(npt,uk,upk, "Heun's Method for N = " + str(N))

    
def rungekutta1(t0,tf,f0,f1,N):
    
    npt = np.linspace(float(t0),float(tf),N)
    
    dt = float((tf-t0)/N)
    uk = np.zeros(N)
    upk = np.zeros(N)
    
    uk[0] = f0
    upk[0] = f1
    
    xK1 = 0.0
    xK2 = 0.0
    
    vK1 = 0.0
    vK2 = 0.0
    
    count = 1
    
    while (count < npt.size):
        xK1 = dt*upk[count-1]
        xK2 = dt*(upk[count-1]+xK1/2)
        vK1 = -dt*uk[count-1]
        vK2 = -dt*(uk[count-1] + vK1/2)
        uk[count] = uk[count-1] + xK2
        upk[count] = upk[count-1] + vK2
        count += 1
        
    plotitokay(npt,uk,upk, "Runge-Kutta 2nd Order for N = " + str(N))
        

def runge2(t0,tf,f0,f1,N):
    
    npt = np.linspace(float(t0),float(tf),N)
    
    dt = float((tf-t0)/N)
    uk = np.zeros(N)
    upk = np.zeros(N)
    
    uk[0] = f0
    upk[0] = f1
    
    xK1 = 0.0
    xK2 = 0.0
    xK3 = 0.0
    xK4 = 0.0
    
    vK1 = 0.0
    vK2 = 0.0
    vK3 = 0.0
    vK4 = 0.0
    
    count = 1
    
    while count < npt.size:
        xK1 = dt*upk[count-1]
        xK2 = dt*(upk[count-1]+xK1/2)
        xK3 = dt*(upk[count-1]+xK2/2)
        xK4 = dt*(upk[count-1]+xK3)
        uk[count] = uk[count-1]+(xK1+2*xK2+2*xK3+xK4)/6
        
        vK1 = -dt*uk[count-1]
        vK2 = -dt*(uk[count-1]+vK1/2)
        vK3 = -dt*(uk[count-1]+vK2/2)
        vK4 = -dt*(uk[count-1]+vK3)
        upk[count] = upk[count-1]+(vK1+2*vK2+2*vK3+vK4)/6
        
        count += 1
        
    plotitokay(npt,uk,upk, "Runge-Kutta 4th Order Method for N = " + str(N))
            
def plotitokay(t,uk,upk):
    fig = plt.figure(figsize = (12,8))
    a = plt.axes()
    
    a.plot(t,uk, 'b', label = "$x(t)$")
    a.plot(t,upk, 'r', label = "$v(t)$")
    
    plt.title(titl)
    a.legend()
    plt.show()