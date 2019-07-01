# ****************************************************************************
# Name           : Main.py
# Author         : Andres Valdez
# Version        : 1.0
# Description    : 3D truss element, pull-back matrix operations
#                  This set of functions receives Nodal displaments, and
#                  provides the transformation matrix, from global to local
# Data	     : 05-05-2019
# ****************************************************************************

from __future__ import unicode_literals
import os, sys
import numpy as np
import scipy as scp
from scipy import optimize
from scipy.linalg import solve
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from math import pi, sin, sqrt, pow

def get_cos(x,y,z,l0):
    '''
    Receives the variations of lengths between two consecutive mechanical
    states. And computes the direction towards the initial state.
    '''
    cx = (x[1] - x[0]) / l0
    cy = (y[1] - y[0]) / l0
    cz = (z[1] - z[0]) / l0
    return cx, cy, cz

def map_global_to_local(cx,cy,cz):
    '''
    Here we define the 3d to 2d mapping. (x,y,z)--->(xo,xf)
    '''
    T = np.zeros((2,6))
    T[0,0], T[0,1], T[0,2] = cx, cy, cz
    T[1,3], T[1,4], T[1,5] = cx, cy, cz
    
    return T

def vec_global_to_local(T,f):
    '''
    Receives any global vector (elemental) and the transformation matrix.
    Returns the local vector for Isoparametric implementation
    '''
    return np.matmul(T,f)
 
def vec_local_to_global(T,f):
    '''
    Receives the local vector and the global to local mapping, and returns
    the global vector.
    '''
    return np.matmul(T.T,f)
   
def gauss_point(npts):
    '''
    This function returns the localization of a gauss points and weights
    consideering the ref domain (-1,1).
    '''
    pg, wg = np.zeros(npts) , np.zeros(npts)
    
    if(npts==1):
        pg[0] = 0.0
        wg[0] = 2
    if(npts==2):
        pg[0] , pg[1] = 1/sqrt(3) , -1/sqrt(3)
        wg[0] , wg[1] = 1 , 1
    return pg, wg

def shape_func_local(npts):
    '''
    Here we define the local shape functions. Isoparametric (-1,1).
    Just for a 2-node line
    shl(0,i,l) = x-derivative of shape function i
    shl(1,i,l) =  shape function defined at node i
    l: gauss point
    '''
    shl = np.zeros((2,2,npts))
    
    pg , wg = gauss_point(npts)
    
    for k in range(npts):
        shl[0,0,k] = -0.5
        shl[1,0,k] = -0.5 * (pg[k] - 1.0)
        
        shl[0,1,k] = 0.5
        shl[1,1,k] = -0.5 * (pg[k] + 1.0)
    
    return shl

def shape_func_global(npts,x):
    '''
    Here we define the global shape functions
    shg(0,i,l) = x-derivative of shape function i
    shg(1,i,l) =  shape function defined at node i
    x[0],x[1] local-coordinates of line nodes
    
    '''
    det = np.zeros(npts)
    shg = np.zeros((2,2,npts))
    
    shl = shape_func_local(npts)
    
    for k in range(npts):
        det[k] = sqrt((x[1]-x[0])**2)
        
        for i in range(2):
            shg[0,i,k] = shl[0,i,k] / det[k]
            shg[1,i,k] = shl[1,i,k]
        
    return shg,det

def pos_truss(shg,x,y,z,l):
    '''
    Here we define the x,y,z-variables for integration @ each gauss-point
    '''
    px = 0.0
    
    for k in range(2):
        px = px + shg[1,k,l] * x[k]
    return px,py,pz

def elemental_force(force,x,npts):
    '''
    Here we define the elemental force vector.
    Receives local forces, local coordinates
    '''
    fe = np.zeros(2)
    pg , wg = gauss_point(npts)
    shg , det = shape_func_global(npts,x)
    for l in range(npts):
        for k in range(2):
            fe[k] = fe[k] + force[k] * wg[l] * shg[1,k,l] * det[l]
        
    return fe

def Green_Lagrange(npts,x,u):
    '''
    Here we define the elemental Green-Lagrange tensor.
    Receives number of gauss points and local coordinates and a local
    displacement.
    '''
    du = np.zeros(2)
    Eu = np.zeros((2,2))
    pg , wg = gauss_point(npts)
    shg , det = shape_func_global(npts,x)
    for l in range(npts):
        for k in range(2):
            du[k] = du[k] + shg[0,k,l] * wg[l] * det[l] * u[k]
    
    Eu[0,0] = du[0] + 0.5 * du[0]**2
    Eu[1,1] = du[1] + 0.5 * du[1]**2
    return Eu

def elemental_int_force(Eu,E,A):
    '''
    Receives the Green-Lagrange strain tensor, the Young modulus, and the
    transversal area. Returns the internal-efforts at both extremal nodes.
    '''
    fi = np.zeros(2)
    fi[0] = E * A * (Eu[0,0] + Eu[0,1])
    fi[1] = E * A * (Eu[1,0] + Eu[1,1])
    #Note that Eu is always diagonal. Its a truss, the kinematics is penalized.
    return fi
    
if __name__ == "__main__":

    #
    # Workflow
    #
    ###################################################################
    ## 0) Definition of a mechanical state
    ###################################################################
    x , y , z = np.array([0,1]) , np.array([0,1]), np.array([0,1]) # This is the initial state
    ###################################################################
    ## 1) From the previous two-solutions, compute the new direction
    ###################################################################
    l = np.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2) # This is the element's length
    cos_theta = get_cos(x,y,z,l)
    print('the angles',type(cos_theta),cos_theta)
    ###################################################################
    ## 2) Evaluate the transformation matrix
    ###################################################################
    A = map_global_to_local(cos_theta[0],cos_theta[1],cos_theta[2])
    print('transformation matrix')
    print(A)
    ###################################################################
    ## 3) Convert the global force in 3d space to a 2d space.
    ###################################################################
    fg = 9.8 * np.array([0,0,-1,0,0,-1]) # This is a hand made external loading
                                         # in our case gravity, negative z direction
    
    fl = vec_global_to_local(A,fg)
    print('global force',fg)
    print('local force',fl)
    ###################################################################
    ## 4) Convert the global coordinates in 3d space to a 2d space.
    ###################################################################
    ug = np.array([x[0],x[1],y[0],y[1],z[0],z[1]])
    ul = vec_global_to_local(A,ug)
    print('global coords',ug)
    print('local coords',ul)
    ###################################################################
    ## 5) Here starts the isoparametric evaluations. Integrate the external-forces
    ###################################################################
    npts = 1 # Number of Gauss points
    int_fl = elemental_force(fl,ul,npts)
    
    print('integrate force',int_fl)
    ###################################################################
    ## 5) Integrate the internal-forces
    ###################################################################
    E , A = 1 , 1
    int_def = Green_Lagrange(npts,ul,2*ul)
    print('The green lagrange strain tensor')
    print(int_def)
    
    int_force = elemental_int_force(int_def,E,A)
    
    print('The internal forces')
    print(int_force)
    
    sys.exit()
    
