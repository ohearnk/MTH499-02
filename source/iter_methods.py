#!/usr/bin/python

from math import fabs
import numpy as np
import numpy.linalg as la

def gauss_seidel(A,b,x0,M,eps):
    if A.shape[0] != A.shape[1]:
        raise InvalidDimensionError

    # x_k and x_{k-1} for computing the stopping
    # criterion ||x_k - x_{k-1}|| < eps
    x_old = np.copy(x0)
    x_new = np.copy(x_old)
    # list of errors as their computed
    e = []
    # iteration number
    k = 0

    while k < M:
        # compute the next iterate
        for i in xrange(A.shape[0]):
            # use the new values as soon as you're computed
            x_new[i] = (b[i] - np.sum(np.delete(A[i,:],i).dot(np.delete(x_new,i)))) / A[i,i]
        # append the error to the list
        e.append(fabs(la.norm(x_new-x_old)))
        # check stopping criterion
        if fabs(la.norm(x_new-x_old)) < eps:
            break
        # save new copy
        x_old = np.copy(x_new)
        # advance
        k += 1
    
    return x_new, e


def jacobi(A,b,x0,M,eps):
    if A.shape[0] != A.shape[1]:
        raise InvalidDimensionError

    # x_k and x_{k-1} for computing the stopping
    # criterion ||x_k - x_{k-1}|| < eps
    x_old = np.copy(x0)
    x_new = np.zeros(b.shape, dtype=float)
    # list of errors as their computed
    e = []
    # iteration number
    k = 0
    
    while k < M:
        # compute the next iterate
        for i in xrange(A.shape[0]):
            x_new[i] = (b[i] - np.sum(np.delete(A[i,:],i).dot(np.delete(x_old,i)))) / A[i,i]
        # append the error to the list
        e.append(fabs(la.norm(x_new-x_old)))
        # check stopping criterion
        if fabs(la.norm(x_new-x_old)) < eps:
            break
        # save new copy
        x_old = np.copy(x_new)
        # advance
        k += 1
    
    return x_new, e


def richardson(A,b,x0,M,eps,omega=1):
    if A.shape[0] != A.shape[1]:
        raise InvalidDimensionError

    # x_k and x_{k-1} for computing the stopping
    # criterion ||x_k - x_{k-1}|| < eps
    x_old = np.copy(x0)
    x_new = np.copy(x_old)
    # residuals
    r = np.zeros(x_old.shape[0])
    # list of errors as their computed
    e = []
    # iteration number
    k = 0

    while k < M:
        # compute the residuals
        for i in xrange(A.shape[0]):
            r[i] = b[i] - np.sum(A[i,:]*x_new)
        # calculate the next iterate 
        x_new = x_old + omega*r
        # append the error to the list
        e.append(fabs(la.norm(x_new-x_old)))
        # check stopping criterion
        if fabs(la.norm(x_new-x_old)) < eps:
            break
        # save new copy
        x_old = np.copy(x_new)
        # advance
        k += 1

    return x_new, e


if __name__ == '__main__':
    # test Richardson method
    A = np.array([[6,1,1],[2,4,0],[1,2,6]])
    b = np.array([12,0,6])
    x0 = np.zeros((3))
    om = 1.0/6.0
#    A = np.array([[1,0.5,0.33],[0.33,1,0.5],[0.5,0.33,1]])
#    b = np.array([0.61,0.61,0.61])
#    x0 = np.zeros((3))

    print 'A ='
    print A
    print 'b =', b
    print 'x_0 =', x0
    print
    
    print "Richardson method:"
    x, e = richardson(A,b,x0,100,0.0001,omega=om)
    print "\tNumber of iterations:", len(e)
    print "\tx_"+str(len(e))+" = ", x
    print "\te = ", e[-1]
    
    # test Jacobi and Gauss-Seidel methods
    A = np.array([[10,-2,3],[1,7,0],[5,2,-11]])
    b = np.array([2,6,-3])
    x0 = np.zeros((3))
    
    print
    print 'A ='
    print A
    print 'b =', b
    print 'x_0 =', x0
    print
    
    print "Jacobi method:"
    x, e = jacobi(A,b,x0,100,0.0001)
    print "\tNumber of iterations:", len(e)
    print "\tx_"+str(len(e))+" = ", x
    print "\te = ", e[-1]

    print
    print "Gauss Seidel method:"
    x, e = gauss_seidel(A,b,x0,100,0.0001)
    print "\tNumber of iterations:", len(e)
    print "\tx_"+str(len(e))+" = ", x
    print "\te = ", e[-1]

