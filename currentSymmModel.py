#
#  migrateBICodesToPy.py
#  
#
#  Created by Claas Hueter on 5/9/12.
#  Copyright (c) 2012 __tinycodes__. All rights reserved.
#

import numpy as np
from scipy.optimize import newton_krylov
from scipy.optimize import anderson
from scipy.optimize import broyden2
from numpy import cosh, zeros_like, mgrid, zeros
import pylab
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad # returns array where at ...[0] the value of the integral is saved
from scipy.integrate import romberg
from scipy.special import *

def integrand_consLaw_LorR(xpVal,xObsVal,yObsVal,Ps_i,Peclet):
	ypVal = Ps_i.__call__(xpVal)
	etaVal = eta(xpVal,xObsVal,ypVal,yObsVal)
	etaSymmVal = eta(-1.*xpVal,xObsVal,ypVal,yObsVal)
	#return 2.*np.exp(-Peclet*( yObsVal - ypVal )) * kv(0,Peclet*etaVal)
	
	if (  abs(xpVal - xObsVal) < 1.e-7 ):
		return 0.
	else:
		return 2.*np.exp(-Peclet*( yObsVal - ypVal )) * (kv(0,Peclet*etaVal)+kv(0,Peclet*etaSymmVal))	
	
def integrate_consLaw__R(xObsVal,yObsVal,Ps_R,Peclet):
	#return romberg(integrand_consLaw_LorR, xProl_rightArray[nAdd], xProl_rightArray[-1], args=(xObsVal,yObsVal,Ps_R,Peclet),rtol=1.48e-04, divmax=10)
	return quad(integrand_consLaw_LorR, xProl_rightArray[nAdd], xProl_rightArray[-1], args=(xObsVal,yObsVal,Ps_R,Peclet), full_output=1)[0]
	
def eta(xp,xObs,yp,yObs):
	return ( (xObs-xp)**2 + (yObs-yp)**2   )**(0.5)

#################################################################################
# parameters
nx_R = 31

nAdd = 5
hx_R = 10./(nx_R-1) 
interpolationFactor = 5
ip = interpolationFactor

kappa = np.zeros(nx_R)
integral = np.zeros(nx_R)

Peclet = 0.2
DeltaVal = ( Peclet * np.pi )**(0.5) * np.exp(Peclet) * erfc( Peclet**(0.5) )
sigma = 1.
angle_phi = 0.15*np.pi   # is the fixed angle by which the ratios of surface tensions at the tip are set 
slope_R = -1./np.tan(angle_phi)

xRightArrayList = [ (ix)*hx_R for ix in range(nx_R)]
xArrayList = xRightArrayList
xSymm_array = np.asarray(xArrayList)

xAddList_R = [ -5.0, -4.0, -3.0, -2.0, -1.0]

xProlRightArrayList = [ -nAdd*hx_R+(ix)*hx_R for ix in range(nx_R+nAdd)]

xProlRightArrayListInterp = [ -nAdd*hx_R+(ix)*hx_R/ip for ix in range((nx_R+nAdd-1)*ip+1)]

xProl_rightArray = np.asarray(xProlRightArrayList)

xProl_rightArrayInterp = np.asarray(xProlRightArrayListInterp)

cntCalls = 0

def residual(P):
	global integral, Peclet, DeltaVal, sigma, angle_phi, slope_R
	global xProl_rightArray, xProl_rightArrayInterp
	global cntCalls
	eqsSystem = np.zeros(nx_R)
	P_R = P
	sigma = P_R[0]
	slope_R = -1./np.tan(angle_phi)
	print 'current guess: ', P_R[0], ' = sigma '
	P_R[0] = 0.
	yAddList_R = ([ el*slope_R for el in xAddList_R ])
	P_Rlist = list(P_R)
	extP_R = np.asarray(yAddList_R+P_Rlist)
	Ps_R = UnivariateSpline(xProl_rightArray, extP_R, s=0)
	P_Rinterp = Ps_R(xProl_rightArrayInterp)
	
	y_p1 = Ps_R.__call__(xProl_rightArray[nAdd+1])
		
	eqsSystem[0] = y_p1-slope_R*xProl_rightArray[nAdd+1]
	
	for iR in range(nx_R-1):
		xObsVal = xProl_rightArray[nAdd+1+iR]
		yObsVal = Ps_R.__call__(xProl_rightArray[nAdd+1+iR])
		yIvVal = -0.5*xObsVal*xObsVal 
		dydxObsVal = Ps_R.__call__(xProl_rightArray[nAdd+1+iR],1)
		d2ydx2ObsVal = Ps_R.__call__(xProl_rightArray[nAdd+1+iR],2)
		dsObsVal = np.sqrt( (1.+dydxObsVal*dydxObsVal) )
		kappaObsVal = d2ydx2ObsVal/(dsObsVal*dsObsVal*dsObsVal) # this way the error on the curvature close to the tip is heavily dependent on the 
																#discretization close to the tip
		integral_R = integrate_consLaw__R(xObsVal, yObsVal,Ps_R,Peclet)
		integral[iR] = Peclet/(2.*np.pi) *  integral_R
		
		#print 'x,y,y_Iv', xObsVal, yObsVal, yIvVal ,'calculated integral', integral[iR], 'analytical value', DeltaVal
		
		eqsSystem[1+iR] = DeltaVal + sigma*kappaObsVal - integral[iR]
	#exit(1)
	return eqsSystem
	
guess = np.asarray([ -0.5*el**2 for el in xArrayList ])
guess[0] = 1.

print 'control proper initial guess: ', guess[0], ' = sigma ' 

sol = newton_krylov(residual, guess, verbose=1, iter=10)
#sol = broyden2(residual, guess, max_rank=50, verbose=1)
#sol = anderson(residual, guess, M=10, verbose=1)
print 'Residual', abs(residual(sol)).max()
#sol = guess

solList = [el for el in list(sol)]
sol_array = np.asarray( solList )

for i in range(len(guess)):
	print solList[i]

for i in range(len(guess)):
	print xSymm_array[i], guess[i]
 
sol_array[0] = 0.

pylab.plot(xSymm_array, sol_array)
pylab.xlabel(' ')
pylab.ylabel(' ')
pylab.title(' ')
pylab.grid(True)
pylab.savefig('testMeiron')
pylab.show()




'''	
def quad_routine(func, a, b, x_list, w_list):
    c_1 = (b-a)/2.0
    c_2 = (b+a)/2.0
    eval_points = map(lambda x: c_1*x+c_2, x_list)
    func_evals = map(func, eval_points)
    return c_1 * sum(array(func_evals) * array(w_list))

def quad_gauss_7(func, a, b):
    x_gauss = [-0.949107912342759, -0.741531185599394, -0.405845151377397, 0, 0.405845151377397, 0.741531185599394, 0.949107912342759]
    w_gauss = array([0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469, 0.381830050505119, 0.279705391489277,0.129484966168870])
    return quad_routine(func,a,b,x_gauss, w_gauss)

def quad_kronrod_15(func, a, b):
    x_kr = [-0.991455371120813,-0.949107912342759, -0.864864423359769, -0.741531185599394, -0.586087235467691,-0.405845151377397, -0.207784955007898, 0.0, 0.207784955007898,0.405845151377397, 0.586087235467691, 0.741531185599394, 0.864864423359769, 0.949107912342759, 0.991455371120813]
    w_kr = [0.022935322010529, 0.063092092629979, 0.104790010322250, 0.140653259715525, 0.169004726639267, 0.190350578064785, 0.204432940075298, 0.209482141084728, 0.204432940075298, 0.190350578064785, 0.169004726639267, 0.140653259715525,  0.104790010322250, 0.063092092629979, 0.022935322010529]
    return quad_routine(func,a,b,x_kr, w_kr)

class Memorize(object):
    def __init__(self, func):
        self.func = func
        self.eval_points = {}
    def __call__(self, *args):
        if args not in self.eval_points:
            self.eval_points[args] = self.func(*args)
        return self.eval_points[args]

def quad(func,a,b):
    ### Output is the 15 point estimate; and the estimated error 
    func = Memorize(func) #  Memoize function to skip repeated function calls.
    g7 = quad_gauss_7(func,a,b)
    k15 = quad_kronrod_15(func,a,b)
    # I don't have much faith in this error estimate taken from wikipedia
    # without incorporating how it should scale with changing limits
    return [k15, (200*scipy.absolute(g7-k15))**1.5]
'''	


##################################################################################
'''
import numpy as np
from scipy.optimize import newton_krylov
from numpy import cosh, zeros_like, mgrid, zeros

# parameters
nx, ny = 75, 75
hx, hy = 1./(nx-1), 1./(ny-1)

P_left, P_right = 0, 0
P_top, P_bottom = 1, 0

def residual(P):
    d2x = zeros_like(P)
    d2y = zeros_like(P)

    d2x[1:-1] = (P[2:]   - 2*P[1:-1] + P[:-2]) / hx/hx
    d2x[0]    = (P[1]    - 2*P[0]    + P_left)/hx/hx
    d2x[-1]   = (P_right - 2*P[-1]   + P[-2])/hx/hx

    d2y[:,1:-1] = (P[:,2:] - 2*P[:,1:-1] + P[:,:-2])/hy/hy
    d2y[:,0]    = (P[:,1]  - 2*P[:,0]    + P_bottom)/hy/hy
    d2y[:,-1]   = (P_top   - 2*P[:,-1]   + P[:,-2])/hy/hy

    return d2x + d2y + 5*cosh(P).mean()**2

# solve
guess = zeros((nx, ny), float)
sol = newton_krylov(residual, guess, verbose=1)
#sol = broyden2(residual, guess, max_rank=50, verbose=1)
#sol = anderson(residual, guess, M=10, verbose=1)
print 'Residual', abs(residual(sol)).max()

# visualize
import matplotlib.pyplot as plt
x, y = mgrid[0:1:(nx*1j), 0:1:(ny*1j)]
plt.pcolor(x, y, sol)
plt.colorbar()
plt.show()
'''
