import math as math
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as sciopt
import random as random

def pythag(x1, y1, z1, x2, y2, z2):
    return math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))

#This removes returns a list of pointmasses brighter than the threshold
def mass_fractalizer(filename, lo_thresh = .01, hi_thresh = 1):
    title = filename[0:-4]+" "+str(100*lo_thresh)+"% to "+str(100*hi_thresh)+"% max brightness"
    
    arr = plt.imread(filename)
    copy = np.copy(arr)
    
    shape = np.shape(arr)
    print("The image has shape", shape)
    
    """
    if(len(shape) == 3):#flatten the array
        new_arr = np.zeros((shape[0], shape[1]))
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(3):
                    new_arr[i, j] = new_arr[i, j] + arr[i, j, k]
        arr = new_arr
    """
        
    min_pix, max_pix = arr[0,0], arr[0,0]
    for i in range(shape[0]):
        for j in range(shape[1]):
            if arr[i,j] < min_pix:
                min_pix = arr[i,j]
            if arr[i,j] > max_pix:
                max_pix = arr[i,j] 
                
    print("The Darkest Pixel is", min_pix)
    print("The Brightest Pixel is", max_pix)
            
    lo_thresh = lo_thresh * max_pix
    hi_thresh = hi_thresh * max_pix
    
    xs, ys, ms = [], [], []
    for i in range(shape[0]):
        for j in range(shape[1]):
            if((arr[i, j] > lo_thresh) and (arr[i, j] <= hi_thresh)):
                xs.append(j)
                ys.append(-i)
                ms.append(arr[i,j])
            else:
                copy[i, j] = 2**16-1
    
    plt.figure()
    plt.imshow(copy, cmap = "gray")
    plt.title(title)
    
    print(len(xs), "pixels survived the thresher")
    return xs, ys, ms

"""
Accepts a list of pointmasses, and returns the ln(C(q,r)) as a list of lists
"""
def Correlation_Integrator(xs, ys, zs, ms, Qs, Rs):
    
    N = len(xs)
    
    #deal with the time
    c = 3.468e-8
    etc = c*len(Rs)*N*N
    hours = int(etc//3600)
    mins = int((etc%3600)//60)
    secs = int(etc%60)
    print("estimated time to completetion is", hours, "hours,", mins, "minutes, and", secs, "seconds")
    print("Please note: the etc is calebrated non-sense.")

    Omega = 0#This is M in the paper, but that got confusing with my Ms variable
    for m in ms:
        Omega += m
    
    Ms = [] #This is an N by R matrix containing the number of points with in a distance r of point p
    for i in range(N):
        Ms.append([])
        for R in Rs:
            Ms[-1].append(ms[i])
    
    for i in range(N):
        for j in range(i+1, N):
            r = pythag(xs[i], ys[i], zs[i], xs[j], ys[j], zs[j])
            for k in range(len(Rs)):
                if(Rs[k] > r):
                    Ms[i][k] = Ms[i][k] + ms[j]
                    Ms[j][k] = Ms[j][k] + ms[i]
    
    lnCqr = [] #the log of the generalized correaltion function as a function of q and r
    
    #now that we have calculated M(p,r), we find Xs and Ys
    for Q in Qs:
        lnCr = [] #the log of the generalized correlation integral as a function of r for some fixed q
        if(Q == 1):#there is a special forumula for Q = 1
            for j, r in enumerate(Rs):
                B = 0
                for i in range(N):
                    M = Ms[i][j]
                    A = math.log(M/Omega)
                    B += ms[i]*A
                lnCr.append(B/Omega)
        else: #Q != 1
            for j, r in enumerate(Rs):
                B = 0
                for i in range(N):
                    M = Ms[i][j]
                    A = math.pow(M,(Q-1))#An intermediate value that makes error propogation simpler
                    B += ms[i]*A#An intermediate value that makes error propogation simpler
                if(B <= 0):
                    print(Q, r)
                lnCr.append((math.log(B) - Q*math.log(Omega))/(Q-1))
        lnCqr.append(lnCr)
    return lnCqr

#Standard linear regression
def lin_reg(xs, ys):
    #plt.figure()
    #plt.scatter(xs, ys)
    Sx = 0
    Sy = 0
    Sxy = 0
    Sxx = 0
    N = 0
    for i in range(len(xs)):
        x = xs[i]
        y = ys[i]
        Sx += x
        Sy += y
        Sxx += x*x
        Sxy += x*y
        N += 1
    #print(N, Sxx, Sx)
    m = (N*Sxy - Sx*Sy)/(N*Sxx-Sx*Sx)
    b = (Sy - m*Sx)/N
    SSres = 0
    SStot = 0
    for i in range(N):
        f = m*x + b
        SSres += (y-f)*(y-f)
        SStot += (y-Sy/N)*(y-Sy/N)
    r2 = 1 - SSres/SStot
    return m, b, r2

def logistic(xs, L, k, a, b):
    return (L-b)/(1+np.exp(k*(xs-a)))+b

def gen_logistic(xs, L, k, a, b, v):
    v = math.exp(v)#this helps with the fitting
    return (L-b)/np.power(1+np.exp(k*(xs-a)), v)+b

#From Falconer 1997
def analytical_multi_sierp_spec(Qs, ps):
    Ds = []
    for Q in Qs:
        if(Q == 1):
            delta = .5
            D = math.pow(ps[0], 1-delta) + math.pow(ps[1], 1-delta) + math.pow(ps[2], 1-delta)
            D = -math.log(D)/math.log(2)/(1-delta-1)
            Ds.append(D)
            D = math.pow(ps[0], 1+delta) + math.pow(ps[1], 1+delta) + math.pow(ps[2], 1+delta)
            D = -math.log(D)/math.log(2)/(1+delta-1)
            Ds[-1] = (Ds[-1] + D)/2
        else:
            D = math.pow(ps[0], Q) + math.pow(ps[1], Q) + math.pow(ps[2], Q)
            D = -math.log(D)/math.log(2)/(Q-1)
            Ds.append(D)
    return Ds

"""
This does curve fitting. 
For xs[i] >= swap, linear regression is prefromed
Otherwise a generalized logistic curve is fitted to the data
And the slope at the inflection point is returned.
"""
def L_sloper(Xs, Yss, Qs, lower, upper, gen = True, swap = 0):
    ms, paramss = [], []
    lower_ind, upper_ind = len(Xs)-1, 0
    for i in range(len(Xs)):
        x = Xs[i]
        if(x > lower and x < upper):
            if i < lower_ind:
                lower_ind = i
            if i > upper_ind:
                upper_ind = i
    for i in range(len(Qs)):
        Ys = Yss[i]
        if(Qs[i] < swap):
            if(gen):
                #First we fit a regular logistic curve
                params, pcov = sciopt.curve_fit(logistic, Xs, Ys, p0 = [-12, 1, 4, -1], bounds = ([-20, 0, 0, -4], [-4, 10, 20, 0]))
                L, k, a, b = params
                #Then we use the parameters of the regular logistic curve as the starting conditions for the general logistic curve
                params, pcov = sciopt.curve_fit(gen_logistic, Xs, Ys, p0 = [L, k, a, b, 0],  bounds = ([-20, 0, 0, -4, -100], [-4, 10, 20, 0, 100]))
                L, k, a, b, v = params
                v = math.exp(v)#This helps with the fitting.
                m = -k*(L-b)*math.pow((v/(v+1)),(v+1))
            else:
                params, pcov = sciopt.curve_fit(logistic, Xs, Ys, p0 = [-12, 1, 4, 0])
                L, k, a, b = params
                m = k*(b-L)/4
            ms.append(m)
            paramss.append(params)
        else:
            #m, b, r2 = lin_reg(Xs, Ys)
            m, b, r2 = lin_reg(Xs[lower_ind:upper_ind], Ys[lower_ind:upper_ind])
            ms.append(m)
            paramss.append((b,))
    return ms, paramss
