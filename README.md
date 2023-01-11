# BsplineTools

***************************************************************************
FILE LIST:
Derivatives of the Bspline basis114.doc  - technical paper
bSplineTools.m - B-spline fitting tools class code
exampleApplication.mlx - Live file demonstrating the use of the class


***************************************************************************
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

......METHODS LIST......

CONSTRUCTOR METHOD:

obj = myBspline(dx,ks,lo,hi);

Arguments:

dx = degree of interpolating polynomial
ks = knot sequence (strictly ascending) in natural units
lo = data low range limit
hi = data hi range limit

Note the B-Spline is only defined on the interval [ lo, hi ]. Allow for
this if it is necessary to extrapolate the fitted function beyond this
range.

CODE METHOD:

Codes data in natural units onto the interval [lo, hi] --> [ 0, 1]. All
fitting and basis function calculation processes are conducted in coded 
units.

xc = obj.code(x)

Arguments:

x = data in natural or engineering units.
 
DECODE METHOD:

Decodes data from the interval [ 0, 1 ] --> [ lo, hi ].

X = obj.code(xc)

Arguments:

xc = data in coded units.

BASIS METHOD:

Calculate basis functions for spline for the current knot sequence.

A = obj.basis(x);   

Arguments:

x = data in natural or engineering units.

DIFFBASIS METHOD:

Calculate the rth derivative of the basis functions with respect to the
input variable

D = obj.diffBasis(x,dr)

Arguments:

x = data in natural or engineering units.
dr = rth derivative to calculate {1}

EVAL METHOD:

Evaluate predictions at values of x. Assumes basis coefficients have been 
identified previously by running the fit method.

Arguments:

x = data in natural or engineering units.

FIT METHOD:

Fit the spline to the data provided. Note obj.a and obj.b must be set prior
to calling this method.

obj = obj.fit(x,y,options,conStructure)

Input Arguments:

x             --> Independent data vector which must lay in the interval [obj.a, obj.b];
y             --> Dependent data.
options       --> optimoptions.fmincon object
conStructure  --> Multi-dimensional structure specifying the
                   constraints
OPTIONS is a optim.options.Fmincon object, that controls the
behaviour of fmincon. To generate this object use options =
optimoptions(@fmincon) at the command line. Then any
user-definable property can be set at the command line using
options.property = value. For example, to display
intermediate results use options.Display = 'iter';

conStructure is a multi-dimensional structure specifiying the
necessary constraints. The structure must have fields:

derivative    --> set to 0,1 or 2 {0} to specify the spline
                  derivative to which the constraint applies.
type          --> set to '==','>=' or '<='
value         --> constraint bound value
x             --> x-ordinates at which constraints apply.
                  Leave empty to specify all training x-ordinates

For example, to specify the constraint that the
minimum prediction from the spline must be >=10 
at x = -2 and x = -1.5, specify:

conStructure.derivative = 0;
conStructure.type = '>=';
conStructure.value = 10;
conStructure.x = [-2;-1.5];

Now assume that a second constraint applies; namely, that the
second derivative of the spline must be negative over all the 
training points. Then set:

conStructure(2).derivative = 2;
conStructure(2).type = '<=';
conStructure(2).value = 0;
conStructure(2).x = [];

CALCDERIVATIVE METHOD:

Calculate the rth derivative of the basis functions at the coordinates
specified.

D = obj.calcDerivative(x,dr)

Arguments:

x = data in natural or engineering units.
dr = rth derivative to calculate {1}

PLOTBASIS METHOD:

A method to plot the B-spline basis functions 

obj.plotBasis

FIND METHOD:

Find occurences of a value over the interval [a,b]

F = obj.find(value, plotflg);

Input Arguments:
 
value     --> Target value
plotflg   --> set to true to plot the results

GENREFKNOTSEQUENCE METHOD:

A method to generate a reference knot sequence for the knot penalty 
function. This sequence is important as it controls the regularisation of
the fitted function. 

obj = obj.genRefKnotSequence(type)

Arguments:

type --> Controls spacing for reference spacing. Set to 0 for linearly 
         spacing, 1 for logarithmic spacing and 2 for reciprical spacing.
         The concept is to use a spacing most consistent with the behaviour
         of the data.

...... PROPERTY LIST ......

d       --> degree of interpolating polynomial
n       --> knot sequence in natural or egineering units
a       --> lower data limit for coding
b       --> upper data limit for coding
alpha   --> spline basis function coefficients
x       --> x-training data
y       --> y-training data
k       --> number of free knots
ak      --> augmented knot sequence in natural units
akc     --> augmented knot sequence in coded units
Gx      --> multiplicative penalty function value
krefc   --> coded reference knot sequence
nb      --> number of basis function
AIC     --> Akaike's Information Criterion
eta     --> constant gain for penalty function (1.1)