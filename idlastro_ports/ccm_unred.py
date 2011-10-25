
from numpy.lib.polynomial import poly1d
        

import numpy

import pylab

def ccm_unred(wave, flux, a_v=None, ebv=None, r_v=3.1):
    """
     NAME:
         CCM_UNRED
     PURPOSE:
         Deredden a flux vector using the CCM 1989 parameterization
     EXPLANATION:
         The reddening curve is that of Cardelli, Clayton, & Mathis (1989 ApJ.
         345, 245), including the update for the near-UV given by O'Donnell
         (1994, ApJ, 422, 158).   Parameterization is valid from the IR to the
         far-UV (3.5 microns to 0.1 microns).
    
         Users might wish to consider using the alternate procedure FM_UNRED
         which uses the extinction curve of Fitzpatrick (1999).
     CALLING SEQUENCE:
         CCM_UNRED, wave, flux, ebv, funred, [ R_V = ]
                 or
         CCM_UNRED, wave, flux, ebv, [ R_V = ]
     INPUT:
         WAVE - wavelength vector (Angstroms)
         FLUX - calibrated flux vector, same number of elements as WAVE
                 If only 3 parameters are supplied, then this vector will
                 updated on output to contain the dereddened flux.
         EBV  - color excess E(B-V), scalar.  If a negative EBV is supplied,
                 then fluxes will be reddened rather than deredenned.
    
     OUTPUT:
         FUNRED - unreddened flux vector, same units & number of elements
                 as FLUX
    
     OPTIONAL INPUT KEYWORD
         R_V - scalar specifying the ratio of total selective extinction
                 R(V) = A(V) / E(B - V).    If not specified, then R_V = 3.1
                 Extreme values of R(V) range from 2.75 to 5.3
    
     EXAMPLE:
         Determine how a flat spectrum (in wavelength) between 1200 A & 3200 A
         is altered by a reddening of E(B-V) = 0.1.   Assume an "average"
         reddening for the diffuse interstellar medium (R(V) = 3.1)
    
           IDL> w = 1200 + findgen(40)*50      ;Create a wavelength vector
           IDL> f = w*0 + 1                    ;Create a "flat" flux vector
           IDL> ccm_unred, w, f, -0.1, fnew  ;Redden (negative E(B-V)) flux vector
           IDL> plot,w,fnew
    
     NOTES:
         (1) The CCM curve shows good agreement with the Savage & Mathis (1979)
                 ultraviolet curve shortward of 1400 A, but is probably
                 preferable between 1200 & 1400 A.
         (2)  Many sightlines with peculiar ultraviolet interstellar extinction
                 can be represented with a CCM curve, if the proper value of
                 R(V) is supplied.
         (3)  Curve is extrapolated between 912 & 1000 A as suggested by
                 Longo et al. (1989, ApJ, 339,474)
         (4) Use the 4 parameter calling sequence if you wish to save the
                   original flux vector.
         (5) Valencic et al. (2004, ApJ, 616, 912) revise the ultraviolet CCM
                 curve (3.3 -- 8.0 um-1).    But since their revised curve does
                 not connect smoothly with longer & shorter wavelengths, it is
                 not included here.
    
     REVISION HISTORY:
           Written   W. Landsman        Hughes/STX   January, 1992
           Extrapolate curve for wavelengths between 900 & 1000 A   Dec. 1993
           Use updated coefficients for near-UV from O'Donnell   Feb 1994
           Allow 3 parameter calling sequence      April 1998
           Converted to IDLV5.0                    April 1998
    """
    # ON_ERROR, 2
    
    
#    if (r_v is None):    
#        r_v = 3.1
    
    x = 10000. / numpy.array(wave)                # Convert to inverse microns 
    npts = x.size
    a = numpy.zeros(npts, dtype=numpy.float)
    b = numpy.zeros(npts, dtype=numpy.float)
    #******************************
    
    #good = numpy.where(ravel(bitwise_and((x > 0.3), (x < 1.1))))[0]       #Infrared
    good = numpy.where( (x >= 0.3) & (x < 1.1) )
    if len(good[0]) > 0:    
        a[good] = 0.574 * x[good] ** (1.61)
        b[good] = -0.527 * x[good] ** (1.61)
    
    #******************************
    
    #good = numpy.where(ravel(bitwise_and((x >= 1.1), (x < 3.3))))[0]           #Optical/NIR
    good = numpy.where( (x >= 1.1) & (x < 3.3) )
    if len(good[0]) > 0:                 #Use new constants from O'Donnell (1994)
        y = x[good] - 1.82
        #     c1 = [ 1. , 0.17699, -0.50447, -0.02427,  0.72085,    $ ;Original
        #                 0.01979, -0.77530,  0.32999 ]               ;coefficients
        #     c2 = [ 0.,  1.41338,  2.28305,  1.07233, -5.38434,    $ ;from CCM89
        #                -0.62251,  5.30260, -2.09002 ]

        #** NOTE **:
        #  IDL poly() wants coefficients starting with A0, then A1 then ...AN where 
        #             AN is the coefficient for X^N
        #             So the coefficients are given in that order
        c1 = numpy.array([1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505])        #from O'Donnell
        c2 = numpy.array([0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347])
        
        #  Numpy's poly1d wants **exactly the opposite order **
        #       so swap 'em

        #stop()
        a[good] = poly1d(c1[::-1])(y)
        b[good] = poly1d(c2[::-1])(y)
    #******************************
    
    good = numpy.where( (x >=3.3) & (x < 8))
    #good = numpy.where(ravel(bitwise_and((x >= 3.3), (x < 8))))[0]           #Mid-UV
    if len(good[0]) > 0:    
        
        y = x[good]
        f_a = numpy.zeros([len(good[0])], dtype=numpy.float)    # f_b = numpy.zeros([ngood], dtype=float32)
        good1 = numpy.where(ravel((y > 5.9)))[0]
        if len(good1[0]) > 0:    
            y1 = y[good1] - 5.9
            f_a[good1] = -0.04473 * y1 ** 2 - 0.009779 * y1 ** 3
            f_b[good1] = 0.2130 * y1 ** 2 + 0.1207 * y1 ** 3
        
        a[good] = 1.752 - 0.316 * y - (0.104 / ((y - 4.67) ** 2 + 0.341)) + f_a
        b[good] = -3.090 + 1.825 * y + (1.206 / ((y - 4.62) ** 2 + 0.263)) + f_b
    
    #   *******************************
    
    #good = numpy.where(ravel(bitwise_and((x >= 8), (x <= 11))))[0]         #Far-UV
    good = numpy.where( (x >= 8) & (x <= 11) )
    if len(good[0]) > 0:    
        y = x[good] - 8.
        c1 = numpy.array([-1.073, -0.628, 0.137, -0.070])
        c2 = numpy.array([13.670, 4.257, -0.420, 0.374])
        a[good] = poly1d(c1[::-1])(y)
        b[good] = poly1d(c2[::-1])(y)
    
    #   *******************************
    #stop()
    
    # Now apply extinction correction to input flux vector
    
    if a_v is None:
        a_v = r_v * ebv

    a_lambda = a_v * (a + b / r_v)
    #print a_v, a, b, r_v, b/r_v
    #print a_lambda
    funred = flux * 10. ** (0.4 * a_lambda)       #Derive unreddened flux
    
    #print "----"
    #print flux
    #print funred
    return funred

