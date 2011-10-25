import numpy
import scipy.interpolate
import scipy.ndimage

""" 
idlcompat.idlbase

Python implementations of various functions from the IDL 
standard library...

This is **VERY INCOMPLETE**  - use at your own risk!

M.D. Perrin, mperrin@ucla.edu
2009-07-01
"""

# TODO better handling for multidimensional arrays here...
# The basic idea is to take the IDL axes spec, pull it into a 
# tuple, then call the Numpy zeros/arange function with the tuple 
# reversed...
def intarr(*l):
    return numpy.zeros(l[::-1],dtype='int')
def fltarr(*l):
    return numpy.zeros(l[::-1],dtype='float32')
def dblarr(*l):
    return numpy.zeros(l[::-1],dtype='float')
def bytarr(*l):
    return numpy.zeros(l[::-1],dtype='byte')
def strarr(*l):
    return numpy.zeros(l[::-1],dtype=numpy.str)


def indgen(*l):
    v = numpy.arange( numpy.asarray(l).prod()  ,dtype='int')
    v.shape=l[::-1]
    return v
def findgen(*l):
    v = numpy.arange( numpy.asarray(l).prod()  ,dtype='float32')
    v.shape=l[::-1]
    return v
def dindgen(*l):
    v = numpy.arange( numpy.asarray(l).prod()  ,dtype='float')
    v.shape=l[::-1]
    return v
def bindgen(*l):
    v = numpy.arange( numpy.asarray(l).prod()  ,dtype='byte')
    v.shape=l[::-1]
    return v


def dist(n):
    y,x = numpy.indices((n,n))
    y-=n/2
    x-=n/2
    r = numpy.sqrt(x**2+y**2)
    return shift(r,n/2,n/2)



def size(array0):

    # try to replicate IDL's behavior with scalars being treated as 
    # 1-element arrays?
    if isinstance(array0,numpy.ndarray): array=array0
    elif isinstance(array0,list): array = numpy.asarray(array0)
    else: array = numpy.asarray([array0])

    # how should this handle axes orders?? 
    # probably best to return in IDL convention x,y,z order

    # TODO convert DTYPE into IDL type codes - 
    # see http://idlastro.gsfc.nasa.gov/idl_html_help/SIZE.html
    # versus NumPy codes from http://idlastro.gsfc.nasa.gov/idl_html_help/SIZE.html
    return [len(array.shape)]+ list(array.shape[::-1]) +['dtype']



def strpos(string, substr,pos=None):
    if pos is not None: return string.find(substr,pos)
    else: return string.find(substr)

def strmid(string, start, length):
    return string[start:start+length]

def strlen(string):
    return len(string)


def plot(x,y,xrange=None,yrange=None,
    xtitle=None,ytitle=None,title=None, xlog=0, ylog=0):
    """ Wrapper for IDL's plot command, including /log axes switches and axes labeling """

    # replicate IDL log axes switches
    if xlog ==0 and ylog==0: pylab.plot(x,y)
    elif xlog==0: pylab.semilogy(x,y)
    elif ylog==0: pylab.semilogx(x,y)
    else: pylab.loglog(x,y)

    #add axis labels and titles
    if xrange is not None: pylab.axis(xmin=xrange[0],xmax=xrange[1])
    if yrange is not None: pylab.axis(ymin=yrange[0],ymay=yrange[1])
    if xtitle is not None: pylab.xlabel(xtitle)
    if ytitle is not None: pylab.ylabel(ytitle)
    if title is not None: pylab.title(title)

def shift(array,s1,s2=0,s3=0):
    dims = array.ndim
    # IDL/Py axes order swappage
    if dims ==1:
        return numpy.roll(array,s1)
    elif dims==2:
        return numpy.roll(numpy.roll(array,s1,axis=1),s2,axis=0)
    elif dims==3:
        return numpy.roll(numpy.roll(numpy.roll(array,s1,axis=2),s2,axis=1),s3,axis=0)
    else:
        print "Don't know how to shift > 3D arrays..."
        return array
    
    
# from http://www.scipy.org/Cookbook/Rebinning, modified slightly
def rebin_avg(a, *newshape0):
    # NOTE: does not swap dimensions for IDL/Py style
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    # swap axes order to allow IDL convention in the newshape0 argument
    newshape = newshape0[::-1]


    shape = a.shape
    lenShape = len(shape)
    factor = numpy.asarray(shape)/numpy.asarray(newshape)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    print ''.join(evList)
    return eval(''.join(evList))

# from http://www.scipy.org/Cookbook/Rebinning, modified slightly
def rebin( a, *newshape0 ):
        '''Rebin an array to a new shape.
        Arbitrary new shape allowed, no interpolation done.
        '''
        # swap axes order to allow IDL convention in the newshape0 argument
        newshape = newshape0[::-1]

        assert len(a.shape) == len(newshape)

        slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
        coordinates = numpy.mgrid[slices]
        indices = coordinates.astype('i')   #choose the biggest smaller integer index
        return a[tuple(indices)]

# from http://www.scipy.org/Cookbook/Rebinning
def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [n.float64, n.float32]:
        a = numpy.cast[float](a)

    m1 = numpy.cast[int](minusone)
    ofs = numpy.cast[int](centre) * 0.5
    old = numpy.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
        return None
    newdims = numpy.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = numpy.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = numpy.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = numpy.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [numpy.arange(i, dtype = numpy.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = numpy.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = numpy.mgrid[nslices]

        newcoords_dims = range(numpy.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (numpy.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print "Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported."
        return None


def bytscl(array, max=None , min=None , nan=0, top=255 ):

    # see http://star.pst.qub.ac.uk/idl/BYTSCL.html
    # note that IDL uses slightly different formulae for bytscaling floats and ints. 
    # here we apply only the FLOAT formula...

    if max is None: max = numpy.nanmax(array)
    if min is None: min = numpy.nanmin(array)
    
    #return (top+0.9999)*(array-min)/(max-min)
    return numpy.maximum(numpy.minimum(  
        ((top+0.9999)*(array-min)/(max-min)).astype(numpy.int16)
        , top),0)

