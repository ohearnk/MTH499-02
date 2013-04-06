#!/usr/bin/python

from math import fabs
import numpy as np
import numpy.linalg as la

def gauss_seidel(A,b,x0,M,eps):
    if A.shape[0] != A.shape[1]:
        raise InvalidDimensionError

    x_old = x0
    x_new = np.zeros(b.shape, dtype=float)
    k = 0

    while k < M:
        for i in xrange(A.shape[0]):
            x_new[i] = (b[i] - np.sum(np.delete(A[i,:],i).dot(np.delete(x_old,i)))) / A[i,i]
        print (k+1), fabs(la.norm(x_new-x_old))
        #print x_old
        #print x_new
        if fabs(la.norm(x_new-x_old)) < eps:
            break
        x_old = np.copy(x_new)
        k += 1

if __name__ == '__main__':
    # test code
    A = np.array([[10,-2,3],[1,7,0],[5,2,-11]])
    b = np.array([2,6,-3])
    x0 = np.zeros((3))

    print 'A =', A
    print 'b =', b
    print 'x_0 =', x0
    print

    gauss_seidel(A,b,x0,100,0.0001)

