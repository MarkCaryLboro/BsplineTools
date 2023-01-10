# BsplineTools

myBspline class written by M. Cary

This class implements the B-spline basis function calculation and fits a  
one-dimensional free knot model to data. Details of the recursive 
calculations for the B-spline basis can be found in the text be DeBoor or
M. Cary's PhD thesis.


******************************VERY IMPORTANT*******************************
The B-spline basis calculations use code from the MBC toolbox. This will
cause difficulties with distribution to clients without an MBC license. 
However, the code can be compiled to mitigate this.
***************************************************************************

Constructor Method:

obj = myBspline(dx,ks,lo,hi);

Arguments:

dx = degree of interpolating polynomial
ks = knot sequence
lo = data low range limit
hi = data hi range limit

Note the B-Spline is only defined on the interval [ lo, hi ]. Allow for
this if it is necessary to extrapolate the fitted function beyond this
range.

code method:

Codes data in natural units onto the interval [lo, hi] --> [ 0, 1]. All
fitting and basis function calculation processes are conducted in coded 
units.

xc = obj.code(x)

Arguments:

x = data in natural or engineering units.
 
decode method:

Decodes data from the interval [ 0, 1 ] --> [ lo, hi ].