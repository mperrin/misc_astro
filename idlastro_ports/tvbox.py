import numpy as np
import matplotlib.pyplot as pl

def tvbox(box=(1,1), xcen=0, ycen=0, center=None,**kwargs):
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
    x = [xcen-box[0], xcen+box[0], xcen+box[0], xcen-box[0], xcen-box[0]]
    y = [ycen-box[1], ycen-box[1], ycen+box[1], ycen+box[1], ycen-box[1]]
    pl.plot(x,y, **kwargs)

