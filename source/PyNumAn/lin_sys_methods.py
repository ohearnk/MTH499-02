#!/usr/bin/python

from math import fabs, sqrt
import numpy as np
import numpy.linalg as la


# iterate using the conjugate directions
def conjugate_gradient(A,b,x0,M,eps,delta):
    if A.shape[0] != A.shape[1]:
        raise InvalidDimensionError

    # initial values
    x = np.copy(x0)
    r = b - np.dot(A,x)
    v = np.copy(r)
    c = np.dot(r,r)
    # list of errors
    e = []
    # iteration number
    k = 0

    while k < M:
        if sqrt(np.dot(v,v)) < delta:
            break
        # save intermediate computations (trade space for compute time)
        z = np.dot(A,v)
        # compute the scalar applied to the descent vector
        t = float(c) / np.dot(v,z)
        # compute the next iterate and residual
        x = x + t*v
        r = r - t*z
        # compute the inner product of r^k
        d = np.dot(r,r)
        # append the error to the list
        e.append(d)
        # check stopping criterion
        if d < eps:
            break
        # update results
        v = r + (float(d)/c)*v
        c = d
        # advance
        k += 1
    return x, e


# fixed number of iterations using the conjugate directions
def conjugate_gradient_direct(A,b,x0,eps):
    if A.shape[0] != A.shape[1]:
        raise InvalidDimensionError

    # inital values
    x = np.copy(x0)
    r = b - np.dot(A,x)
    v = np.copy(r)
    # list of errors
    e = []
    # iteration number
    k = 0

    while k < A.shape[0]:
        # compute new values
        u = np.dot(A,v)
        t = float(np.dot(v,r)) / np.dot(v,u)
        x = x + t*v
        r = b - np.dot(A,x)
        # compute the norm of r^k
        d = la.norm(r)
        # append the error to the list
        e.append(d)
        # stop if residual norm is sufficiently small
        if d < eps:
            break
        # update results
        s = -np.dot(r,u) / np.dot(v,u)
        v = r + s*v
        # advance
        k += 1
    
    return x, e


# solve the system Ax = b by equivalently solving
# the quadratic optimization problem:
#   min ||Ax - b||^2
# for the quadratic form
#   Q(x) = ||Ax - b||^2 = x^TAx - 2x^Tb = <x,Ax>-2<x,b>
def gradient_descent(A,b,x0,M,eps):
    if A.shape[0] != A.shape[1]:
        raise InvalidDimensionError

    # inital guess
    x = np.copy(x0)
    # list of errors
    e = []
    # iteration number
    k = 0

    while k < M:
        # compute the residual
        r = b - A.dot(x)
        # compute the scalar applied to the descent vector
        t = float(r.dot(r)) / r.dot(A.dot(r))
        # compute the next iterate
        x = x + t*r
        # stop if the residual norm is sufficiently small
        d = la.norm(r)
        # append the error to the list
        e.append(r)
        # check stopping criterion
        if d < eps:
            break
        # advance
        k += 1
    
    return x, e


def gauss_seidel(A,b,x0,M,eps):
    if A.shape[0] != A.shape[1]:
        raise InvalidDimensionError

    # x_k and x_{k-1} for computing the stopping
    # criterion ||x_k - x_{k-1}|| < eps
    x_old = np.copy(x0)
    x_new = np.copy(x_old)
    # list of errors
    e = []
    # iteration number
    k = 0

    while k < M:
        # compute the next iterate
        for i in xrange(A.shape[0]):
            # use the new values as soon as you're computed
            x_new[i] = (b[i] - np.sum(np.delete(A[i,:],i).dot(np.delete(x_new,i)))) / float(A[i,i])
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
    # list of errors
    e = []
    # iteration number
    k = 0
    
    while k < M:
        # compute the next iterate
        for i in xrange(A.shape[0]):
            x_new[i] = (b[i] - np.sum(np.delete(A[i,:],i).dot(np.delete(x_old,i)))) / float(A[i,i])
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
    # list of errors
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
