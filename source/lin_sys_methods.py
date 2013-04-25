#!/usr/bin/python

from math import fabs
import numpy as np
import numpy.linalg as la

### Auxilary methods ###

# real-valued inner-product
def inner_product(x,y):
    return np.sum(x.dot(y))


def is_symmetric(A):
    for i in xrange(A.shape[0]):
        for j in range(i+1,A.shape[0]):
            if A[i,j] != A[j,i]:
                return False
    return True

# TODO: implement Cholesky factorization for 
# checking if a matrix is positive definite


### Methods for solving linear systems ###

# solve the system Ax = b by equivalently solving
# the quadratic optimization problem:
#   min ||Ax - b||^2
# for the quadratic form
#   Q(x) = ||Ax - b||^2 = x^TAx - 2x^Tb = <x,Ax>-2<x,b>
def gradient_descent(A,b,x0,M,eps):
    if A.shape[0] != A.shape[1]:
        raise InvalidDimensionError

    # x^k and x^{k-1} for computing the stopping
    # criterion ||x^k - x^{k-1}|| < eps
    x_old = np.copy(x0)
    x_new = np.copy(x_old)
    # list of errors as their computed
    e = []
    # iteration number
    k = 0

    while k < M:
        # compute the residual
        r = b - A.dot(x_old)
        # compute the scalar applied to the descent vector
        t = inner_product(r,r) / inner_product(r,A.dot(r))
        # compute the next iterate
        x_new = x_old + t*r
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
