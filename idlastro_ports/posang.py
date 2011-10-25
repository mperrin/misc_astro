import numpy as np


def posang(ra1, dc1, ra2, dc2, units_in='radians', units_out='radians'):
    """
    Ported from IDL posang.pro by Marshall Perrin

     NAME:
           POSANG
     PURPOSE:
           Computes rigorous position angle of source 2 relative to source 1

     EXPLANATION:
           Computes the rigorous position angle of source 2 (with given RA, Dec)
           using source 1 (with given RA, Dec) as the center.

     CALLING SEQUENCE:
           POSANG, U, RA1, DC1, RA2, DC2, ANGLE

     INPUTS:
           U    -- Describes units of inputs and output:
                   0:  everything radians
                   1:  RAx in decimal hours, DCx in decimal
                           degrees, ANGLE in degrees
           RA1  -- Right ascension of point 1
           DC1  -- Declination of point 1
           RA2  -- Right ascension of point 2
           DC2  -- Declination of point 2

       OUTPUTS:
           ANGLE-- Angle of the great circle containing [ra2, dc2] from
                   the meridian containing [ra1, dc1], in the sense north
                   through east rotating about [ra1, dc1].  See U above
                   for units.

       PROCEDURE:
           The "four-parts formula" from spherical trig (p. 12 of Smart's
           Spherical Astronomy or p. 12 of Green' Spherical Astronomy).

       EXAMPLE:
           For the star 56 Per, the Hipparcos catalog gives a position of
           RA = 66.15593384, Dec = 33.94988843 for component A, and
           RA = 66.15646079, Dec =  33.96100069 for component B.   What is the
           position angle of B relative to A?

           IDL> RA1 = 66.15593384/15.d   & DC1 = 33.95988843
           IDL> RA2 = 66.15646079/15.d   & DC2 = 33.96100069
           IDL> posang,1,ra1,dc1,ra2,dc2, ang
                will give the answer of ang = 21.4 degrees
       NOTES:
           (1) If RA1,DC1 are scalars, and RA2,DC2 are vectors, then ANGLE is a
           vector giving the position angle between each element of RA2,DC2 and
           RA1,DC1.   Similarly, if RA1,DC1 are vectors, and RA2, DC2 are scalars,
           then DIS is a vector giving the position angle of each element of RA1,
           DC1 and RA2, DC2.    If both RA1,DC1 and RA2,DC2 are vectors then ANGLE
           is a vector giving the position angle between each element of RA1,DC1
           and the corresponding element of RA2,DC2.    If then vectors are not the
           same length, then excess elements of the longer one will be ignored.

           (2) Note that POSANG is not commutative -- the position angle between
            A and B is theta, then the position angle between B and A is 180+theta
       PROCEDURE CALLS:
            ISARRAY()
       HISTORY:
           Modified from GCIRC, R. S. Hill, RSTX, 1 Apr. 1998
           Use V6.0 notation W.L. Mar 2011

    """

    scalar_ = not (hasattr(ra1, '__len__') or hasattr(ra2, '__len__'))
    if scalar_:
        if (ra1 == ra2) and  (dc1 == dc2):
            angle = 0.0e0
            print 'Positions are equal:  ', ra1, dc1
            return angle

    d2r = np.pi / 180.0
    h2r = np.pi / 12.0

    
    if units_in=='radians':
        rarad1 = ra1
        rarad2 = ra2
        dcrad1 = dc1
        dcrad2 = dc2
    elif units_in=='hours':
        rarad1 = ra1 * h2r
        rarad2 = ra2 * h2r
        dcrad1 = dc1 * d2r
        dcrad2 = dc2 * d2r
    elif units_in=='degrees':
        rarad1 = ra1 * d2r
        rarad2 = ra2 * d2r
        dcrad1 = dc1 * d2r
        dcrad2 = dc2 * d2r
    else:
        raise ValueError("units_in must be one of {radians, hours, degrees} to define the units for RA.")
    

    radif = rarad2 - rarad1
    angle = np.arctan2(np.sin(radif), np.cos(dcrad1) * np.tan(dcrad2) - np.sin(dcrad1) * np.cos(radif))


    if units_out =='degrees':
        angle = angle / d2r
    elif units_out == 'radians':
        pass
    else:
        raise ValueError('units_out must be one of {radians, degrees}')
 
    return angle


