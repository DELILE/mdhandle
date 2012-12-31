# mdhandle,  http://github.com/dtlussier/mdhandle
# Copyright (c) 2008, Dan Lussier @ Oxford University, FBG Group
# Released under the GNU General Public License, v2
"""
Credit: Method 2 in http://www.scipy.org/Cookbook/Least_Squares_Circle using LSQ with jacobian.

TODO: Document CircleFit and test_fit() as other parts of mdhandle

""" 

import numpy as np
from scipy import optimize

class CircleFit(object):
    """
    Fits input data with circle using LSQ method with Jacobian.
    """

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def calc_R(self, xc, yc):
        """ 
        Calculate the distance of each 2D points from the centre c=(xc, yc)
    
        """
        return np.sqrt((self.x-xc)**2 + (self.y-yc)**2)


    def f(self, c):
        """
        Calculate the algebraic distance between the 2D points and the mean
        circle centred at c=(xc, yc).
    
        """
        Ri = self.calc_R(*c)
    
        return Ri - Ri.mean()

    def Df(self, c):
        """ 
        Jacobian of f.
    
        The axis corresponding to derivatives must be coherent with the 
        col_deriv option of leastsq.
    
        """
        xc, yc     = c
        df2b_dc    = np.empty((len(c), self.x.size))

        Ri = self.calc_R(xc, yc)
        df2b_dc[ 0] = (xc - self.x)/Ri                   # dR/dxc
        df2b_dc[ 1] = (yc - self.y)/Ri                   # dR/dyc
        df2b_dc     = df2b_dc - df2b_dc.mean(axis=1)[:, np.newaxis]

        return df2b_dc

    def fit_circle(self, basename='circle'):

        # coordinates of the barycentre
        x_m = np.mean(self.x)
        y_m = np.mean(self.y)

        centre_estimate = x_m, y_m
        self.centre, self.ier = optimize.leastsq(self.f, centre_estimate, 
                                                 Dfun=self.Df, col_deriv=True)

        xc, yc = self.centre
        Ri        = self.calc_R(xc, yc)
        self.R         = Ri.mean()
        self.residu    = sum((Ri - self.R)**2)
        self.residu2   = sum((Ri**2-self.R**2)**2)


# -----------------------------------------------------------------------------

x_test=np.array([  9, 35, -13,  10,  23,   0])
y_test = np.array([ 34, 10,   6, -14,  27, -10])

def test_fit(x=x_test, y=y_test):
    """ 
        Test and plot fitting

    """
    import matplotlib.pyplot as plt
    
    fitting = CircleFit(x, y)
    fitting.fit_circle()

    xc, yc = fitting.centre

    f = plt.figure( facecolor='white')  #figsize=(7, 5.4), dpi=72,
    plt.axis('equal')

    theta_fit = np.linspace(-np.pi, np.pi, 180)

    x_fit = xc + fitting.R*np.cos(theta_fit)
    y_fit = yc + fitting.R*np.sin(theta_fit)
    plt.plot(x_fit, y_fit, 'b-', lw=2)

    plt.plot([xc], [yc], 'bD', mec='y', mew=1)

    # draw
    plt.xlabel('x')
    plt.ylabel('y')

    plt.plot(x, y, 'ro', label='data', ms=8, mec='b', mew=1)
    plt.draw()
