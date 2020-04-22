#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:31:44 2020

@author: nickvalverde
"""

import numpy as np
from numpy import cos, sin, pi
import matplotlib.pyplot as plt

#--Initial settings
mm = 1e-3
mr = 1e-3
emittance = 1.*mm*mm
Np = 100 #Number of particles
sigma_0 = pi/3 #phase angle. This angle gives a f = 4 to match single particle thin lens code.
d = 2. #[m] #full drift length
Lp = 2*d #[m] lattice period

#--Calculate focal length f through stability criterion
f = d*np.sqrt(2/(1-cos(sigma_0))) #focal length

#--Functions to calcluate Twiss parameters
def alpha_func(d,f, sigma_0):
    numerator = -d/f
    denominator = np.sqrt(sin(sigma_0))
    return numerator/denominator

def beta_func(d, f, sigma_0):
    numerator = 2*d - d**3/(4*f**2)
    denominator = sin(sigma_0)
    return numerator/denominator

def gamma_func(d, f, sigma_0):
    numerator = sin(sigma_0) + d**2/f**2
    denominator = 2*d - d**3/(4*f**2)
    return numerator/denominator



def ellipse_load(d, f, sigma_0, emittance, Np, fill = False, match = True):
    #--Initialze x and x' arrays to be populated
    x = []
    xprime = []
    count = 0 #Counter for while loops

    #--Check Match condition
    if match == True:
        #--Calculate twiss parameters along with w and w'
        alpha, betta, gamma = alpha_func(d, f, sigma_0), beta_func(d, f, sigma_0), gamma_func(d, f, sigma_0)
        w = np.sqrt((2.*d - d**3/(4*f**2))/sin(sigma_0))
        wp = -(d/f)/np.sqrt(2*d - d**3/(4*f**2))
        #--Use derived w and w' to evaluate x and x'
        if fill == False: # Perimeter load
             while count < Np:
                theta = 2*pi*count/Np
                xpos = np.sqrt(emittance)*w*cos(theta)
                xprime_pos = np.sqrt(emittance)*sin(theta)/w + np.sqrt(emittance)*wp*cos(theta)

                x.append(xpos)
                xprime.append(xprime_pos)
                count += 1
        else: # Uniform load
            while count < Np:
             #--Fill  Procedure
                up = np.random.random()
                upsi = np.random.random()
                psi = 2*pi*upsi
                xpos = np.sqrt(emittance*np.sqrt(up))*w*cos(psi)
                xprime_pos = np.sqrt(emittance*np.sqrt(up))*sin(psi)/w + np.sqrt(emittance*np.sqrt(up))*wp*cos(psi)

                x.append(xpos)
                xprime.append(xprime_pos)
                count +=1
        #--Output twiss parameters
        print("alpha = %f" %alpha)
        print("betta = %f[m]" %betta)
        print("gamma = %f" %gamma)


    #--If matched conditions is not desired use random number generator
        #to be coded later.
    else:
        pass

    #--Turn x and xprime into np arrays
    x = np.array(x)
    xprime = np.array(xprime)


    return x, xprime

#-Create Transport matrix M
M0 = np.array( [[1.,d/2.],[0,1.]] )         # drift matrix
Mf = np.array( [[1.,0],[-1./f,1.]] )     # focus matrix
Md = np.array( [[1.,0],[ 1./f,1.]] )     # defocus matrix
#M = M0 @ Md @ M0 @ M0 @ Mf @ M0
M = M0 @ Mf @ M0 @ M0 @ Md @ M0
def advance(x, xprime, Matrix):
    x_new = []
    xprime_new = []

    #Loop through x and xprime and advance particles. Append advanced coordinates.
    for i in range(len(x)):
        a = x[i]
        b = xprime[i]
        vec = np.array([[a], [b]])
        vec = Matrix@vec
        x_new.append(vec[0][0])
        xprime_new.append(vec[1][0])

    x, xprime = np.array(x_new), np.array(xprime_new)
    return x, xprime


x, xprime = ellipse_load(d, f, sigma_0, emittance, Np, fill = False, match = True)

fig, ax = plt.subplots(figsize=(7,7))
ax.scatter(x/mm, xprime/mm, s = 0.5)
ax.set_title(r"d = %g[m], f = %g[m],$\sigma_0 = $%g[rad] " %(d, f, sigma_0), fontsize = 18)
ax.set_xlabel(r'$x$[mm]', fontsize = 18)
ax.set_ylabel(r"x'[mr]", fontsize = 18)
xmax = max(x)/mm
xpmax = max(xprime)/mr
x_axis, yx = [0,0],[0,xmax]
xp_axis, ypx = [0,0], [0,xpmax]
ax.set_title(r"d = %g[m], f = %g[m],$\sigma_0 = $%g[rad] " %(d, f, sigma_0), fontsize = 18)
ax.plot(yx, x_axis, c = 'k', lw = .75)
ax.plot(xp_axis, ypx, c = 'k', lw = .75)
x_annote = '(%1.3f, %g)' %(xmax, 0)
xp_annote = '(%g, %1.3f)' %(0, xpmax)
ax.annotate(x_annote, (xmax, 0), xytext = (xmax-xmax/5,xpmax/20 ))
ax.annotate(xp_annote, (0,xpmax), xytext = (xmax/35, xpmax))


count = 0

while count < 200:
    x, xprime = advance(x,xprime,M)
    count += 1

fig, ax = plt.subplots(figsize=(7,7))
ax.scatter(x/mm, xprime/mm, s = 0.5)
ax.set_title(r"Pass 51, d = %g[m], f = %g[m],$\sigma_0 = $%g[rad] " %(d, f, sigma_0), fontsize = 18)
ax.set_xlabel(r'$x$[mm]', fontsize = 18)
ax.set_ylabel(r"x'[mr]", fontsize = 18)
xmax = max(x)/mm
xpmax = max(xprime)/mr
x_axis, yx = [0,0],[0,xmax]
xp_axis, ypx = [0,0], [0,xpmax]
ax.set_title(r"d = %g[m], f = %g[m],$\sigma_0 = $%g[rad] " %(d, f, sigma_0), fontsize = 18)
ax.plot(yx, x_axis, c = 'k', lw = .75)
ax.plot(xp_axis, ypx, c = 'k', lw = .75)
x_annote = '(%1.3f, %g)' %(xmax, 0)
xp_annote = '(%g, %1.3f)' %(0, xpmax)
ax.annotate(x_annote, (xmax, 0), xytext = (xmax-xmax/5,xpmax/20 ))
ax.annotate(xp_annote, (0,xpmax), xytext = (xmax/35, xpmax))