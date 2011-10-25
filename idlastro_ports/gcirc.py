import numpy as np

def gcirc(ra1, dc1, ra2, dc2, units_in='degrees', units_out='arcsec'):
    """
    Ported from IDL gcirc.pro by Marshall Perrin
    
     NAME:
         GCIRC
     PURPOSE:
         Computes rigorous great circle arc distances.
     EXPLANATION:
         Input position can either be either radians, sexigesimal RA, Dec or
         degrees.   All computations are double precision.
    
     CALLING SEQUENCE:
          GCIRC, U, RA1, DC1, RA2, DC2, DIS
    
     INPUTS:
          U    -- integer = 0,1, or 2: Describes units of inputs and output:
                  0:  everything radians
                  1:  RAx in decimal hours, DCx in decimal
                           degrees, DIS in arc seconds
                  2:  RAx and DCx in degrees, DIS in arc seconds
          RA1  -- Right ascension or longitude of point 1
          DC1  -- Declination or latitude of point 1
          RA2  -- Right ascension or longitude of point 2
          DC2  -- Declination or latitude of point 2
    
     OUTPUTS:
          DIS  -- Angular distance on the sky between points 1 and 2
                  See U above for units;  double precision
    
     PROCEDURE:
          "Haversine formula" see
          http://en.wikipedia.org/wiki/Great-circle_distance
    
     NOTES:
           (1) If RA1,DC1 are scalars, and RA2,DC2 are vectors, then DIS is a
           vector giving the distance of each element of RA2,DC2 to RA1,DC1.
           Similarly, if RA1,DC1 are vectors, and RA2, DC2 are scalars, then DIS
           is a vector giving the distance of each element of RA1, DC1 to
           RA2, DC2.    If both RA1,DC1 and RA2,DC2 are vectors then DIS is a
           vector giving the distance of each element of RA1,DC1 to the
           corresponding element of RA2,DC2.    If the input vectors are not the
           same length, then excess elements of the longer ones will be ignored.
    
           (2) The function SPHDIST provides an alternate method of computing
            a spherical distance.
    
           (3) The haversine formula can give rounding errors for antipodal
           points.
    
     PROCEDURE CALLS:
          None
    
       MODIFICATION HISTORY:
          Written in Fortran by R. Hill -- SASC Technologies -- January 3, 1986
          Translated from FORTRAN to IDL, RSH, STX, 2/6/87
          Vector arguments allowed    W. Landsman    April 1989
          Prints result if last argument not given.  RSH, RSTX, 3 Apr. 1998
          Remove ISARRAY(), V5.1 version        W. Landsman   August 2000
          Added option U=2                      W. Landsman   October 2006
          Use double precision for U=0 as advertised R. McMahon/W.L.  April 2007
          Use havesine formula, which has less roundoff error in the
                 milliarcsecond regime      W.L. Mar 2009
    """

    
    d2r = np.pi / 180.0
    as2r = np.pi / 648000
    h2r = np.pi / 12.
    
    # Convert input to double precision radians
    if units_in == 'radians':
        rarad1 = numpy.asarray(ra1, dtype=float64)
        rarad2 = numpy.asarray(ra2, dtype=float64)
        dcrad1 = numpy.asarray(dc1, dtype=float64)
        dcrad2 = numpy.asarray(dc2, dtype=float64)
    elif units_in == 'hours':    
        rarad1 = ra1 * h2r
        rarad2 = ra2 * h2r
        dcrad1 = dc1 * d2r
        dcrad2 = dc2 * d2r
    elif units_in == 'degrees':    
        rarad1 = ra1 * d2r
        rarad2 = ra2 * d2r
        dcrad1 = dc1 * d2r
        dcrad2 = dc2 * d2r
    else:    
        raise ValueError("units_in must be one of {radians, hours, degrees} to define the units for RA.")
    
    
    deldec2 = (dcrad2 - dcrad1) / 2.0
    delra2 = (rarad2 - rarad1) / 2.0
    sindis = np.sqrt(np.sin(deldec2) * np.sin(deldec2) + np.cos(dcrad1) * np.cos(dcrad2) * np.sin(delra2) * np.sin(delra2))
    dis = 2.0 * np.arcsin(sindis)
    
    if units_out =='arcsec':
        dis = dis / as2r
    elif units_out =='degrees':
        dis = dis / d2r
    elif units_out == 'radians':
        pass
    else:
        raise ValueError('units_out must be one of {radians, degrees, arcsec}')
    
    return dis


