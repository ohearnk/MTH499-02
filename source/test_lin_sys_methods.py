from matplotlib import pyplot
import numpy as np

import PyNumAn.lin_sys_methods as lsm

def main():
    # test Richardson method
    A = np.array([[1,0.5,0.33],[0.33,1,0.5],[0.5,0.33,1]])
    b = np.array([0.61,0.61,0.61])
    x0 = np.zeros((3))

    print 'A ='
    print A
    print 'b =', b
    print 'x_0 =', x0
    print
    
    print "Richardson method:"
    x, e = lsm.richardson(A,b,x0,100,0.0001)
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
    x, e = lsm.jacobi(A,b,x0,100,0.0001)
    print "\tNumber of iterations:", len(e)
    print "\tx_"+str(len(e))+" = ", x
    print "\te = ", e[-1]

    print
    print "Gauss Seidel method:"
    x, e = lsm.gauss_seidel(A,b,x0,100,0.0001)
    print "\tNumber of iterations:", len(e)
    print "\tx_"+str(len(e))+" = ", x
    print "\te = ", e[-1]

    # test gradient descent and conjugate gradient methods
    A = np.array([[4,-1,1],[-1,4,-2],[1,-2,4]])
    b = np.array([12,-1,5])
    x0 = np.array([0,0,0])

    print
    print 'A ='
    print A
    print 'b =', b
    print 'x_0 =', x0
    print

    print "Gradient descent method:"
    x, e = lsm.gradient_descent(A, b, x0, 200, 1e-4)
    print "\tNumber of iterations:", len(e)
    print "\tx_"+str(len(e))+" = ", x
    print "\te = ", e[-1]
    
    print "Direct Conjugate gradient method:"
    x, e = lsm.conjugate_gradient_direct(A, b, x0, 1e-4)
    print "\tNumber of iterations:", len(e)
    print "\tx_"+str(len(e))+" = ", x
    print "\te = ", e[-1]

    print "Conjugate gradient method:"
    x, e = lsm.conjugate_gradient(A, b, x0, 200, 1e-4, 1e-4)
    print "\tNumber of iterations:", len(e)
    print "\tx_"+str(len(e))+" = ", x
    print "\te = ", e[-1]


if __name__ == '__main__':
    main()
