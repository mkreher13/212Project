# 22.212 Project: Nodal Expansion Method with CMFD
# Miriam Kreher (2018)

# Imports
import numpy as np
from diff_opts import *
from material import *
from cross_sections import *
from construct import *
from solver import *
from plotter import *
import copy

# Variables
options = DiffusionOpts1D()
options.read("input1.inp")
nBins = options.numBins
nGrps = options.numGroups
nu = 2.4
A = np.zeros((nBins, 4))
ERR = []
BoundERR = []
AERR = []

E = [0,1]
f1 = [0,0]
df1 = [0,0]
f2 = [0,0]
df2 = [0,0]
f3 = [0,0]
df3 = [0,0]
f4 = [0,0]
df4 = [0,0]

f1[0] = 2.*E[0]-1.
df1[0] = 2.
f2[0] = 6.*E[0]*(1.-E[0])-1.
df2[0] = 6.-12.*E[0]
f3[0] = 6.*E[0]*(1.-E[0])*(2.*E[0]-1.)
df3[0] = -6.*(6.*E[0]**2.-6.*E[0]+1.)
f4[0] = 6.*E[0]*(1.-E[0])*(5.*E[0]**2.-5.*E[0]+1.)
df4[0] = -120.*E[0]**3.+180.*E[0]**2.-72.*E[0]+6.

f1[1] = 2.*E[1]-1.
df1[1] = 2.
f2[1] = 6.*E[1]*(1.-E[1])-1.
df2[1] = 6.-12.*E[1]
f3[1] = 6.*E[1]*(1.-E[1])*(2.*E[1]-1.)
df3[1] = -6.*(6.*E[1]**2.-6.*E[1]+1.)
f4[1] = 6.*E[1]*(1.-E[1])*(5.*E[1]**2.-5.*E[1]+1.)
df4[1] = -120.*E[1]**3.+180.*E[1]**2.-72.*E[1]+6.

M = Material()
M.read()

XS = CrossSections(options)
XS.problem1(options, M.data)
Matrix=Construct(options)
# print(XS.D[0])
# print(XS.removal[0])
# print(XS.fis[0])
# print(XS.abs[0])

sol = Solve(options)

# Iterations
for j in range(0,10):

    ####################################
    # Solving CMFD
    # Where A is the matrix of diffusion coefficients
    # F is the source matrix
    Matrix.constructA(options, XS.D[0], XS.removal[0], XS.Gscat[0], nu*XS.fis[0], j)
    # print(Matrix.A)
    # print(Matrix.F)
    Matrix.invert(Matrix.A)
    sol.solve(options, Matrix.inv, Matrix.F)
    # print(sol.flux)
    # print(Matrix.Ds)
    if j == 0:
        FirstFlux = copy.copy(sol.flux)
    LastFlux = copy.copy(sol.flux)
    
    ####################################
    # Solving NEM
    # Where C is the coefficient matrix,
    # BC is the boundary condition version of C,
    # and S is the RHS of those equations
    for n in range(0,nBins):

        if n == 0:
        	# Solving for both BC at once in this function
        	Matrix.constructBC(options, n, sol.flux, sol.k, XS.D[0], XS.removal[0], nu*XS.fis[0], 
                E, f1, df1, f2, df2, f3, df3, f4, df4)
        else:
        	Matrix.constructC(options, n, sol.flux, sol.k, XS.D[0], XS.removal[0], nu*XS.fis[0], 
                E, f1, df1, f2, df2, f3, df3, f4, df4)
        Matrix.invert(Matrix.C)
        sol.coefs(Matrix.inv, Matrix.S)
        
        sol.current(n, sol.a, options, XS.D[0], XS.removal[0], sol.flux, Matrix,
            E, f1, df1, f2, df2, f3, df3, f4, df4)
        # print(sol.a)
        A[n,:] = sol.a
    sol.node_bal_check(options.delta, sol.flux, XS.removal[0], nu*XS.fis[0], A, Matrix.Ds, Matrix.Dh, sol.k)

    Matrix.Dh[0] = -sol.J[0]/sol.flux[0] 
    Matrix.Dh[len(sol.J)-1] = -sol.J[len(sol.J)-1]/sol.flux[len(sol.flux)-1] 
    # print(sol.J[:])
    # print(Matrix.Ds[0]*sol.flux[0]+sol.J[0], -Matrix.Ds[len(sol.J)-1]*sol.flux[len(sol.flux)-1]+sol.J[len(sol.J)-1])
    
    ####################################
    # Plotting results
    ERR.append(max(sol.ERR))
    BoundERR.append(max(sol.BoundERR))
    AERR.append(max(sol.AERR))

    results = Plotter()
results.plot_flux(options, sol.flux, j)
results.plot_flux_change(options, FirstFlux, LastFlux)
results.plot_node_err(ERR, BoundERR, AERR)
    # print (FirstFlux[:]-LastFlux[:])/LastFlux[:]









