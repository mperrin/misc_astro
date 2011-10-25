import numpy as np
import matplotlib.pyplot as pl


def tvcircle(radius=1, xcen=0, ycen=0, center=None,**kwargs):
    """
        draw a circle on an image.

            radius
            xcen
            ycen
            center=     tuple in (Y,X) order.
    """
    if center is not None:
        xcen=center[1]
        ycen=center[0]
    t = np.arange(0, pi * 2.0, 0.01).reshape((len(t), 1))
    x = radius * np.cos(t) + xcen
    y = radius * np.sin(t) + ycen
    pl.plot(x,y, **kwargs)

