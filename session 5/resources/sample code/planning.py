import numpy as np
from math import pi,asin,sin
import lattice_utils as lu
#import matplotlib.pyplot as plt
from  mpl_toolkits.axes_grid.grid_helper_curvelinear import GridHelperCurveLinear
from mpl_toolkits.axes_grid.axislines import Subplot

'''import matplotlib.pyplot as plt

def dynamic_range(Efixed,E,E_max,theta_range = [10,120],step = 10, color = 'k',showplot = True):

    #modify to allow fixed Ef or fixed Ei, and input of scattering
    #angles
    omega = np.linspace(0,E_max,100)
    theta_s = np.arange(theta_range[0]*np.pi/180,theta_range[1]*np.pi/180,step*np.pi/180)
    Q = np.empty([theta_s.size,omega.size],float)

    if Efixed == "Ef": 
        kf = np.sqrt((E)/2.072)
        ki = np.sqrt((omega+E)/2.072)
    elif Efixed =="Ei":
        ki = np.sqrt((E)/2.072)
        kf = np.sqrt((E-omega)/2.072)

    for i, theta in enumerate(theta_s):

        Q[i] = np.sqrt(ki**2 + kf**2 - 2*ki*kf*np.cos(theta))
        
        if showplot:
            plt.plot(Q[i],omega,lw=1,ls = '--',color = color)
            txt = "$2\\theta_s$ = {0}$^o$".format(np.round(theta*180/np.pi,1))
            plt.text(Q[i,i+25],omega[50],txt, bbox=dict(fc = '1',lw = 0,alpha = 0.05),rotation = 75,color = color)
            plt.xlabel('Q ($\\AA^{-1}$)')
            plt.ylabel('Energy Transfer (meV)')
            title = 'Accessible dynamic range for {1} fixed = {0} meV'.format(E,Efixed)
            plt.title(title)
            plt.grid(True)
            plt.show()

    return omega,Q

'''
def spec_twoTheta(Efixed,E,E_T,Q):
    
    if Efixed == "Ef": 
        kf = np.sqrt((E)/2.072)
        ki = np.sqrt((E_T+E)/2.072)
    elif Efixed =="Ei":
        ki = np.sqrt((E)/2.072)
        kf = np.sqrt((E-E_T)/2.072)

    theta =np.arccos(-(Q**2 - ki**2 - kf**2)/ki/kf/2.)
    
    return theta*180./np.pi


def Bragg_angle(wavelength,q,rlatt):

    d = lu.dspacing(q,rlatt) 
    print(wavelength/d/2)
    tth = 360./pi*asin(wavelength/d/2)
    
    print('\n')
    print( '\t wavelength = {:.2f}, Q = [{:.2f} {:.2f} {:.2f}], Two-theta = {:.3f}'.format(wavelength,q[0],q[1],q[2],tth))
    return

def TOF_par(q,tth,rlatt):

    d = lu.dspacing(q,rlatt)
    wavelength = 2*d*sin(tth*pi/360)
    E = (9.044/wavelength)**2
    k = 2*pi/wavelength
    velocity = 629.62*k # m/s
    print('Q = [{:.2f} {:.2f} {:.2f}]\n d = {:.3f} \n Two-theta = {:.2f}\n wavelength = {:.3f} Angstrom\n Energy = {:3f} meV\n Velocity = {:3f} m/s'.format(q[0],q[1],q[2],d,tth,wavelength,E,velocity)) 
    
    return d,wavelength,E,velocity

def Recip_space(sample):

    """
    Set up general reciprocal space grid for plotting Miller indicies in a general space.
    Would be cool if returned fig object had a custom transformation so that all data added
    to plot after it has been created can be given in miller indicies
    grid for custom transform.
    """
    def tr(x, y):
        x, y = np.asarray(x), np.asarray(y)
        return x, y-x

    def inv_tr(x,y):
        x, y = np.asarray(x), np.asarray(y)
        return x, y+x

    grid_helper = GridHelperCurveLinear((tr, inv_tr))    
    
    fig = plt.figure(1, figsize=(7, 4))
    ax = Subplot(fig, 1, 1, 1, grid_helper=grid_helper)

    rlatt = sample.star_lattice
    [xs,ys,zs] = sample.StandardSystem

    fig.add_subplot(ax)


    ax.grid(True)
    return


def Al_peaks(wavelength = 1.0):
    
    energy = (9.044/wavelength)**2
    # Al_lattice 
    a = 4.0498; b = 4.0498; c = 4.0498
    aa =90; bb = 90; cc = 90
    
    latt = lu.lattice(a,b,c,aa,bb,cc)  
    rlatt = lu.recip_lattice(latt)
    # Al is FCC so peaks must be all even or all odd
    peaks = [[1,1,1],[2,0,0],[2,2,0],[3,1,1],[2,2,2],[4,0,0],[3,3,1],[4,2,0],
            [4,2,2],[5,1,1],[3,3,3],[4,4,0],[5,3,1],[4,4,2],[6,0,0]]
    
    print('\tAluminum Bragg Peaks')
    print('\tNeutron wavelength {:.2f} Angstroms ({:.2f} meV)\n'.format(wavelength,energy))
    print('\t H   K   L\tQ (AA-1)    d(AA)     2theta     2theta(l/2)     2theta(l/3)')
    print('\t----------------------------------------------------------------------------')
    
    for p in peaks:
        modQ = lu.modVec(p,rlatt)
        dsp = lu.dspacing(p,rlatt)
        
        if abs(wavelength/(2*dsp)) < 1:
            tth = 360/pi * asin(wavelength/(2*dsp))
            line = '\t {:d}   {:d}   {:d}\t {:.3f}      {:.3f}      {:.2f}\t   {:.2f}\t   {:.2f}'
        else :
            tth = 'NaN'
            line = '\t {:d}   {:d}   {:d}\t {:.3f}      {:.3f}      {:s}\t   {:.2f}\t   {:.2f}'
        
        if abs((wavelength/2)/(2*dsp)) < 1:
            tth2 = 360/pi * asin((wavelength/2)/(2*dsp))
        else :
            tth2 = 'NaN'
            line = '\t {:d}   {:d}   {:d}\t {:.3f}      {:.3f}      {:s}\t   {:s}\t\t   {:.2f}'
        
        if abs((wavelength/3)/(2*dsp)) < 1:
            tth3 = 360/pi * asin((wavelength/3)/(2*dsp))
        else :
            tth3 = 'NaN'
            line = '\t {:d}   {:d}   {:d}\t {:.3f}      {:.3f}      {:s}\t   {:s}\t\t   {:s}'
        # else :
            # tth = 360/pi * asin(wavelength/(2*dsp))
            # tth2 = 360/pi * asin((wavelength/2)/(2*dsp))
            # tth3 = 360/pi * asin((wavelength/3)/(2*dsp))
        
        # line = '\t {:d}   {:d}   {:d}\t {:.3f}      {:.3f}      {:.2f}\t   {:.2f}\t   {:.2f}'
        print(line.format(p[0],p[1],p[2],modQ,dsp,tth,tth2,tth3))

    return


if __name__ == '__main__':

    Al_peaks()