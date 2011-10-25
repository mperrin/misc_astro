import numpy as np



def sunpos(jd=None, return_all=False, radian=False):
    """
    sunpos: Compute the RA and Dec of the sun on a given date. 
    Converted from IDL astro sunpos.pro to Python by Marshall Perrin. 


    doc string of IDL sunpos follows: 


     NAME:
           SUNPOS
     PURPOSE:
           To compute the RA and Dec of the Sun at a given date.
    
     CALLING SEQUENCE:
           SUNPOS, jd, ra, dec, [elong, obliquity, /RADIAN ]
     INPUTS:
           jd    - The Julian date of the day (and time), scalar or vector
                   usually double precision
     OUTPUTS:
           ra    - The right ascension of the sun at that date in DEGREES
                   double precision, same number of elements as jd
           dec   - The declination of the sun at that date in DEGREES
    
     OPTIONAL OUTPUTS:
           elong - Ecliptic longitude of the sun at that date in DEGREES.
           obliquity - the obliquity of the ecliptic, in DEGREES
    
     OPTIONAL INPUT KEYWORD:
           /RADIAN - If this keyword is set and non-zero, then all output variables
                   are given in Radians rather than Degrees
    
     NOTES:
           Patrick Wallace (Rutherford Appleton Laboratory, UK) has tested the
           accuracy of a C adaptation of the sunpos.pro code and found the
           following results.   From 1900-2100 SUNPOS  gave 7.3 arcsec maximum
           error, 2.6 arcsec RMS.  Over the shorter interval 1950-2050 the figures
           were 6.4 arcsec max, 2.2 arcsec RMS.
    
           The returned RA and Dec are in the given date's equinox.
    
           Procedure was extensively revised in May 1996, and the new calling
           sequence is incompatible with the old one.
     METHOD:
           Uses a truncated version of Newcomb's Sun.    Adapted from the IDL
           routine SUN_POS by CD Pike, which was adapted from a FORTRAN routine
           by B. Emerson (RGO).
     EXAMPLE:
           (1) Find the apparent RA and Dec of the Sun on May 1, 1982
    
           IDL> jdcnv, 1982, 5, 1,0 ,jd      ;Find Julian date jd = 2445090.5
           IDL> sunpos, jd, ra, dec
           IDL> print,adstring(ra,dec,2)
                    02 31 32.61  +14 54 34.9
    
           The Astronomical Almanac gives 02 31 32.58 +14 54 34.9 so the error
                   in SUNPOS for this case is < 0.5".
    
           (2) Find the apparent RA and Dec of the Sun for every day in 1997
    
           IDL> jdcnv, 1997,1,1,0, jd                ;Julian date on Jan 1, 1997
           IDL> sunpos, jd+ dindgen(365), ra, dec    ;RA and Dec for each day
    
     MODIFICATION HISTORY:
           Written by Michael R. Greason, STX, 28 October 1988.
           Accept vector arguments, W. Landsman     April,1989
           Eliminated negative right ascensions.  MRG, Hughes STX, 6 May 1992.
           Rewritten using the 1993 Almanac.  Keywords added.  MRG, HSTX,
                   10 February 1994.
           Major rewrite, improved accuracy, always return values in degrees
           W. Landsman  May, 1996
           Added /RADIAN keyword,    W. Landsman       August, 1997
           Converted to IDL V5.0   W. Landsman   September 1997
    """

    if jd is None:
        import datetime
        dd = datetime.datetime(2011,1,1,0)
        now = dd.utcnow()
        jd = gd2jd(now.year, now.month, now.day, now.hour+(now.minute+now.second/60.)/60.)
    jd = np.asarray(jd)

    
    dtor = np.pi / 180.0       #(degrees to radian, double precision)
    
    #  form time in Julian centuries from 1900.0
    
    t = (jd - 2415020.0e0) / 36525.0e0
    
    #  form sun's mean longitude
    
    l = (279.696678 + np.mod((36000.768925 * t) , 360.0)) * 3600.0
    
    #  allow for ellipticity of the orbit (equation of centre)
    #  using the Earth's mean anomaly ME
    
    me = 358.475844 + ((35999.049750 * t) % 360.0e0)
    ellcor = (6910.1 - 17.2 * t) * np.sin(me * dtor) + 72.3 * np.sin(2.0 * me * dtor)
    l = l + ellcor
    
    # allow for the Venus perturbations using the mean anomaly of Venus MV
    
    mv = 212.603219 + np.mod((58517.803875 * t) , 360.0)
    vencorr = 4.8 * np.cos((299.1017 + mv - me) * dtor) + 5.5 * np.cos((148.3133 + 2.0 * mv - 2.0 * me) * dtor) + 2.5 * np.cos((315.9433 + 2.0 * mv - 3.0 * me) * dtor) + 1.6 * np.cos((345.2533 + 3.0 * mv - 4.0 * me) * dtor) + 1.0 * np.cos((318.15 + 3.0 * mv - 5.0 * me) * dtor)
    l = l + vencorr
    
    #  Allow for the Mars perturbations using the mean anomaly of Mars MM
    
    mm = 319.529425 + ((19139.858500 * t) % 360.0e0)
    marscorr = 2.0 * np.cos((343.8883 - 2.0 * mm + 2.0 * me) * dtor) + 1.8 * np.cos((200.4017 - 2.0 * mm + me) * dtor)
    l = l + marscorr
    
    # Allow for the Jupiter perturbations using the mean anomaly of
    # Jupiter MJ
    
    mj = 225.328328 + ((3034.6920239 * t) % 360.0e0)
    jupcorr = 7.2 * np.cos((179.5317 - mj + me) * dtor) + 2.6 * np.cos((263.2167 - mj) * dtor) + 2.7 * np.cos((87.1450 - 2.0 * mj + 2.0 * me) * dtor) + 1.6 * np.cos((109.4933 - 2.0 * mj + me) * dtor)
    l = l + jupcorr
    
    # Allow for the Moons perturbations using the mean elongation of
    # the Moon from the Sun D
    
    d = 350.7376814 + ((445267.11422 * t) % 360.0e0)
    mooncorr = 6.5 * np.sin(d * dtor)
    l = l + mooncorr
    
    # Allow for long period terms
    
    longterm = +6.4 * np.sin((231.19 + 20.20 * t) * dtor)
    l = l + longterm
    l = (l + 2592000.0e0) % 1296000.0e0
    longmed = l / 3600.0e0
    
    # Allow for Aberration
    
    l = l - 20.5e0
    
    # Allow for Nutation using the longitude of the Moons mean node OMEGA
    
    omega = 259.183275 - np.mod((1934.142008 * t) , 360.0)
    l = l - 17.2 * np.sin(omega * dtor)
    
    # Form the True Obliquity
    
    oblt = 23.452294 - 0.0130125 * t + (9.2 * np.cos(omega * dtor)) / 3600.0e0
    
    # Form Right Ascension and Declination
    
    l = l / 3600.0e0
    #print l, oblt
    ra = np.arctan2(np.sin(l * dtor) * np.cos(oblt * dtor), np.cos(l * dtor))

    if hasattr(ra, '__len__'):
        # ra is an array
        neg = np.where(ra < 0.0e0)
        ra[neg] = ra[neg] + 2.0 * np.pi
    else:
        # ra is a scalar
        if ra < 0: ra += 2.0*np.pi
        
    dec = np.arcsin(np.sin(l * dtor) * np.sin(oblt * dtor))
    
    if radian :    
        oblt = oblt * dtor
        longmed = longmed * dtor

    else:
        ra = ra / dtor
        dec = dec / dtor
   
    if return_all: 
        return (ra, dec, oblt, longmed)
    else: 
        return (ra, dec)
        




