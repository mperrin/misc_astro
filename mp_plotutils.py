import numpy as np
import pylab as p

import matplotlib


# purely cosmetic: print e.g. '100' instead of '10^2' for axis labels
class NicerLogFormatter(matplotlib.ticker.LogFormatter):
    """ A LogFormatter subclass to print nicer axis labels
        e.g. '100' instead of '10^2'

        Parameters
        ----------
        threshhold : int
            Absolute value of log base 10 above which values are printed in exponential notation.
            By default this is 3, which means you'll get labels like 10^-4, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10^4 ...

        usage:
          ax = gca()
          ax.set_yscale("log")
          ax.set_xscale("log")
          ax.xaxis.set_major_formatter(NiceLogFormatter())
    """
    def __init__(self, threshhold=3):
        self.threshhold = threshhold
    def __call__(self,val,pos=None):
        if abs(np.log10(val)) > self.threshhold:
            return "$10^{%d}$" % np.log10(val)
        elif val >= 1:
            return "%d"%val
        else:
            return "%.1f"%val

def savepdf(filename, **kwargs):
    p.savefig(filename, transparent=True, **kwargs)
    os.system("open "+filename)

