# Purpose: This program employs several iterative
#   methods for finding the roots of single-variate,
#   real-valued functions.  The methods available
#   for use include:
#       + bisection method
#       + Newton's method
#       + Steffenson's method
#       + secant method

import argparse
import csv
from itertools import izip_longest
from math import cos, sin, tan, pi, exp, fabs, sqrt
from texttable import Texttable


def newtons_method(M, eps, delta, f, fp, x_0):
    # record the values of x_n, f(x_n) at each iteration 
    x_n = [ x_0 ]
    f_n = [ f(x_0) ]

    if fabs(f_n[0]) < eps:
        return x_n, f_n
        
    for i in xrange(1,M):
        x_n.append( x_n[i-1] - f_n[i-1]/fp(x_n[i-1]) )
        f_n.append( f(x_n[i]) )

        if fabs(x_n[i]-x_n[i-1]) < delta or fabs(f_n[i]) < eps:
            break
        
    return x_n, f_n


def horners_method(f_coeffs, x_i):
    f_i = f_coeffs[-1]
    fp_i = 0.0

    for i in range(len(f_coeffs)-1,-1,-1):
        fp_i = f_i+x_i*fp_i
        f_i = f_coeffs[i] + x_i*f_i
    
    return f_i, fp_i

def newtons_method_hm(M, eps, delta, f_coeffs, x_0):
    # record the values of x_n, f(x_n) at each iteration 
    x_n = [ x_0 ]
    f_0, fp_0 = horners_method(f_coeffs, x_0)
    f_n = [ f_0 ]

    if fabs(f_n[0]) < eps:
        return x_n, f_n
        
    for i in xrange(1,M):
        f_i, fp_i = horners_method(f_coeffs, x_n[i-1])
        x_n.append( x_n[i-1] - f_i/fp_i )
        f_n.append( f_i )

        if fabs(x_n[i]-x_n[i-1]) < delta or fabs(f_n[i]) < eps:
            break
        
    return x_n, f_n


# create command-line argument parser object
parser = argparse.ArgumentParser(description="Execute the prescribed iterative methods")

# add accepted positional arguments
parser.add_argument("methods",
                    nargs='*', type=str, default='newton',
                    choices=['newton','newton_hm'],
                    help="Iterative methods to run")

# add accepted optional arguments
parser.add_argument("-c", "--csv",
                    metavar='CSVFILE', nargs=1, type=str,
                    default=None,
                    help="Output iterations to CSV file.")
parser.add_argument("-f", "--func_vals", action='store_true',
                    help="Show function values for the iterates of each methods")
parser.add_argument("-v", "--verbosity", action='count',
                    help="Increase output verbosity")

# parse arguments
args = parser.parse_args()

# method parameters
f = lambda x: 3.0*x**10 - 72.5*x**8 - 24.3*x**7 - 5.2*x**4 - 5.2*pi*x - 14.2
fp = lambda x: 30.0*x**9 - 580.0*x**7 - 170.1*x**6 - 20.8*x**3 - 5.2*pi
f_coeffs = [ 3.0, 0.0, -72.5, -24.3, 0.0, 0.0, -5.2, 0.0, 0.0, -5.2*pi, -14.2 ]
M = 100
x_0 = 1.0
eps = 0.0001
delta = 0.0001

# list of lists of iterates produced by the methods
iterates = [ ]

# texttable attributes
header = [ 'Iteration' ]

if 'newton' in args.methods:
    x_n, f_n = newtons_method(M,eps,delta,f,fp,x_0)

    if args.verbosity >= 1:
        print "\nNewton's method:"
        print "\tNumber of iterations: " + str(len(x_n))
        print "\tx_"+str(len(x_n)-1)+" = "+str(x_n[-1])
        print "\tf(x_"+str(len(f_n)-1)+") = "+str(f_n[-1])
    
    iterates.append(x_n)
    header.append("Newton's: x_n")
    
    if args.func_vals:    
        iterates.append(f_n)
        header.append("Newton's: f(x_n)")

if 'newton_hm' in args.methods:
    x_n, f_n = newtons_method_hm(M,eps,delta,f_coeffs,x_0)

    if args.verbosity >= 1:
        print "\nNewton's method (Horner):"
        print "\tNumber of iterations: " + str(len(x_n))
        print "\tx_"+str(len(x_n)-1)+" = "+str(x_n[-1])
        print "\tf(x_"+str(len(f_n)-1)+") = "+str(f_n[-1])
    
    iterates.append(x_n)
    header.append("Newton's (Horner): x_n")
    
    if args.func_vals:    
        iterates.append(f_n)
        header.append("Newton's (Horner): f(x_n)")

# create a texttable and add the records
table = Texttable()
# set the table style
table.set_deco(Texttable.HEADER)
# set precision of floating point data type
table.set_precision(7)
# set column data types in table
table.set_cols_dtype(['i']+['f' for i in iterates])
# set the table column alignment
table.set_cols_align(['r']+['r' for i in iterates])
# add table column headers
table.header(header)

# add the records to the table
i = 0
for ROW in izip_longest(*iterates):
    table.add_row([i]+list(ROW))
    i=i+1

# draw table
print table.draw()

if args.csv != None:
    # create CSV writer
    with open(args.csv, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',
            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    #    writer.writerow(header)
    #    for RECORD in records:
    #        writer.writerow(RECORD)

