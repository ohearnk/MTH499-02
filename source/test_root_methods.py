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

import PyNumAn.root_methods as rm

def main():
    # create command-line argument parser object
    parser = argparse.ArgumentParser(description="Execute the prescribed iterative methods")

    # add accepted positional arguments
    parser.add_argument("methods",
                        nargs='*', type=str, default='bisection',
                        choices=['bisection','halley','newton','secant','secant_swap','steffenson'],
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
    f = lambda x: x**2-20.0
    fp = lambda x: 2.0*x
    M = 1000
    a = 3
    b = 3.5
    x_0 = 3.0
    x_1 = 5.0
    eps = 0.0001
    delta = 0.0001

    #f = lambda x: 0.0005*exp(x) + tan(x) - cos(x**2-2.0)
    #fp = lambda x: 0.0005*exp(x) + 1.0/(cos(sqrt(x))**2) + 2.0*x*sin(x**2-2.0)
    #M = 1000
    #a = 3
    #b = 3.5
    #x_0 = 3.0
    #x_1 = 5.0
    #eps = 0.0001
    #delta = 0.0001

    #f = lambda x: (4.4*x**5 - exp(cos(x - 5.0)))/sqrt(5.0*pi*x)
    #fp = lambda x: (0.126257*exp(cos(5.0-x)-0.252313*x*sin(5.0-x)*exp(cos(5.0-x))+4.9958*x**5)/(x**(3.0/2.0)))
    #fpp = lambda x: 2.0
    #M = 1000
    #a = 0.1
    #b = 5.0
    #x_0 = 0.1
    #x_1 = 1.0
    #eps = 0.0001
    #delta = 0.0001

    #f = lambda x: (x*(x*(x*(x*(x+200)+275)+312)+380)+400)-1750
    #M = 1000
    #x_0 = 0.0
    #x_1 = 1.0
    #eps = 0.0001
    #delta = 0.0001

    # list of lists of iterates produced by the methods
    iterates = [ ]

    # texttable attributes
    header = [ 'Iteration' ]

    if 'bisection' in args.methods:
        a_n, b_n = rm.bisection_method(M,eps,delta,f,a,b)

        if args.verbosity >= 1:
            print "Bisection method:"
            print "\tNumber of iterations: " + str(len(a_n))
            print "\ta_"+str(len(a_n)-1)+"= "+str(a_n[-1])
            print "\tb_"+str(len(b_n)-1)+"= "+str(b_n[-1])
            print "\tf(a_"+str(len(a_n)-1)+") = "+str(f(a_n[-1]))
            print "\tf(b_"+str(len(a_n)-1)+") = "+str(f(b_n[-1]))

        iterates.append(a_n)
        iterates.append(b_n)
        header.append('Bisection: a_n')
        header.append('Bisection: b_n')
        
        if args.func_vals:
            # TODO: fix method to store lists of function values
            pass 
    #        iterates.append(fa_n)
    #        iterates.append(fb_n)
    #        header.append('Bisection: f(a_n)')
    #        header.append('Bisection: f(b_n)')

    if 'newton' in args.methods:
        x_n, f_n = rm.newtons_method(M,eps,delta,f,fp,x_0)

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

    if 'halley' in args.methods:
        x_n, f_n = rm.halleys_method(M,eps,delta,f,fp,fpp,x_0)

        if args.verbosity >= 1:
            print "\nHalley's method:"
            print "\tNumber of iterations: " + str(len(x_n))
            print "\tx_"+str(len(x_n)-1)+" = "+str(x_n[-1])
            print "\tf(x_"+str(len(f_n)-1)+") = "+str(f_n[-1])
        
        iterates.append(x_n)
        header.append("Halley's: x_n")
        
        if args.func_vals:    
            iterates.append(f_n)
            header.append("Halley's: f(x_n)")

    if 'steffenson' in args.methods:
        x_n, f_n = rm.steffensons_method(M,eps,delta,f,x_0)

        if args.verbosity >= 1:
            print "\nSteffenson's method:"
            print "\tNumber of iterations: " + str(len(x_n))
            print "\tx_"+str(len(x_n)-1)+" = "+str(x_n[-1])
            print "\tf(x_"+str(len(f_n)-1)+") = "+str(f_n[-1])
        
        iterates.append(x_n)
        header.append("Steffenson's: x_n")
        
        if args.func_vals:    
            iterates.append(f_n)
            header.append("Steffenson's: f(x_n)")

    if 'secant' in args.methods:
        x_n, f_n = rm.secant_method(M,eps,delta,f,x_0,x_1)

        if args.verbosity >= 1:
            print "\nSecant method:"
            print "\tNumber of iterations: " + str(len(x_n))
            print "\tx_"+str(len(x_n)-1)+" = "+str(x_n[-1])
            print "\tf(x_"+str(len(f_n)-1)+") = "+str(f_n[-1])
        
        iterates.append(x_n)
        header.append('Secant: x_n')
        
        if args.func_vals:    
            iterates.append(f_n)
            header.append('Secant: f(x_n)')
        
    if 'secant_swap' in args.methods:
        x_n, f_n = rm.secant_method_swap(M,eps,delta,f,x_0,x_1)

        if args.verbosity >= 1:
            print "\nSecant method (swap):"
            print "\tNumber of iterations: " + str(len(x_n))
            print "\tx_"+str(len(x_n)-1)+" = "+str(x_n[-1])
            print "\tf(x_"+str(len(f_n)-1)+") = "+str(f_n[-1])
        
        iterates.append(x_n)
        header.append('Secant swap: x_n')

        if args.func_vals:    
            iterates.append(f_n)
            header.append('Secant swap: f(x_n)')

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
            writer.writerow(header)
            for RECORD in records:
                writer.writerow(RECORD)


if __name__ == "__main__":
    main()
