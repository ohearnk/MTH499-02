# Purpose: This programs implements the bisection
#   method of iteratively finding the roots of functions
#   in a given interval [a,b]

import math

def bisection_method(M, eps, delta, f, a, b):
    # record a_n and b_n
    a_n = [ a ]
    b_n = [ b ]
    # number of iterations
    m = 0

    # no root in [a,b], return immediately
    if f(a)*f(b) > 0.0:
        print "No root in ["+str(a)+",  "+str(b)+"]."
        return a_n, b_n, m

    # iterate
    u = f(a)
    v = f(b)
    e = b - a
    c = a + (b - a)/2.0
    w = f(c)

    while u*v < 0.0 and m < M \
    and math.fabs(e) >= delta \
    and math.fabs(w) >= eps:
        # [a,c] contains root
        if u*w < 0.0:
            # update b and v = f(b)
            b = c
            v = w
        # [c,b] contains root
        else:
            # update a and u = f(a)
            a = c
            u = w
        # update for next iteration
        e = b - a
        c = a + (b - a)/2.0
        w = f(c)
        # save the iterations
        a_n.append(a)
        b_n.append(b)
        # increase the number of iterations
        m = m + 1
    
    return a_n, b_n, m


# main program
f = lambda x: x**2-2.0
M = 100
a = 1
b = 6
eps = 0.000000001
delta = 0.000000001

a_n, b_n, m = bisection_method(M,eps,delta,f,a,b)

print "a_n:",a_n[-1]
print "b_n:", b_n[-1]
print "f(a_"+str(m)+") = "+str(f(a_n[-1]))
print "f(b_"+str(m)+") = "+str(f(b_n[-1]))
print "\nNumber of iterations: " + str(m)
