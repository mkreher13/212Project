#class to construct linear system of equations

import numpy as np

class Construct():
	
###############################################################
	
	def __init__(self, opt):
		self.Dh = np.zeros(opt.numBins+1)
		self.Ds = np.zeros(opt.numBins+1)
		
###############################################################
		
	def constructA(self, options, D, rXS, Gscat, XSf, j):

		nGrps = options.numGroups
		nBins = options.numBins
		rank = nBins*nGrps
		delta = options.delta

		self.A = np.zeros((rank, rank))
		self.F = np.zeros((rank, rank))

		for g in range(0, nGrps):
			for x in range(0,nBins):
				i = g*nBins+x

				# Calculate the average diffusion theory between 
				# current cell and the adjacent cell
				# and account for vacuum boundary conditions
				if x == 0:
					DAdj = 1
					deltaAdj = 0 # 0 for Zero flux, 4 for Marshak
				else:
					DAdj = D
					deltaAdj = delta
				dLeft = (2*D*DAdj)/(deltaAdj*D+delta*DAdj)
				DhLeft = self.Dh[i]
				
				if x == (nBins-1):
					DAdj = 1
					deltaAdj = 0 # 0 for Zero flux, 4 for Marshak
				else:
					DAdj = D
					deltaAdj = delta
				dRight = (2*D*DAdj)/(deltaAdj*D+delta*DAdj)
				DhRight = self.Dh[i+1]

				self.Ds[i] = dLeft
				if x == (nBins-1):
					self.Ds[i+1] = dRight

				self.A[i,i] = dLeft+dRight+rXS*delta-DhRight+DhLeft
				# print(self.A)
				# print dLeft+dRight

				# if j !=0:
				# 	if x == 0:
				# 		# print j, x
				# 		self.A[i,i] = dRight+rXS*delta-DhRight+DhLeft
				# 	elif x == (nBins-1):
				# 		# print j, x
				# 		self.A[i,i] = dLeft+rXS*delta-DhRight+DhLeft

				if x != 0:
					self.A[i,i-1] = -dLeft+DhLeft
				if x != (nBins-1):
					self.A[i,i+1] = -dRight-DhRight

				if g == 0:
					self.F[i,i] = XSf*delta

		# print(self.A)
		# print(self.F)
				
###############################################################
	def invert(self, A):
		
		self.inv=np.linalg.inv(A)

###############################################################
		
	def constructC(self, opt, n, flux, k, D, XSr, XSf, 
		E, f1, df1, f2, df2, f3, df3, f4, df4):

		XSd = D/opt.delta**2
		a3 = -5.*XSr/(3.*(60.*XSd+XSr))
		a4 = 7*XSr/(420.*(XSd+3.*XSr))

		self.C = np.zeros((4, 4))
		self.C[0,:] = [f1[1]+f3[1]*a3, f2[1]+f4[1]*a4, -f1[0]-a3*f3[0], -f2[0]-a4*f4[0]]
		self.C[1,:] = [df1[1]+df3[1]*a3, df2[1]+df4[1]*a4, -df1[0]-df3[0]*a3, -df2[0]-df4[0]*a4]
		self.C[2,:] = [-df1[1]-df3[1]*a3+df1[0]+df3[0]*a3, -df2[1]-df4[1]*a4+df2[0]+df4[0]*a4, 0, 0]
		self.C[3,:] = [0, 0, -df1[1]-df3[1]*a3+df1[0]+df3[0]*a3, -df2[1]-df4[1]*a4+df2[0]+df4[0]*a4]

		self.S = np.zeros(4)
		self.S[0] = flux[n]-flux[n-1]
		self.S[1] = 0.
		self.S[2] = (XSf/k-XSr)*opt.delta**2*flux[n-1]/D
		self.S[3] = (XSf/k-XSr)*opt.delta**2*flux[n]/D

###############################################################

	def constructBC(self, opt, n, flux, k, D, XSr, XSf, 
		E, f1, df1, f2, df2, f3, df3, f4, df4):

		XSd = D/opt.delta**2
		self.C = np.zeros((4, 4))
		self.S = np.zeros(4)

		a3 = -5.*XSr/(3.*(60.*XSd+XSr))
		a4 = 7*XSr/(420.*(XSd+3.*XSr))

		# Vacuum boundary condition
		# self.C[0,:] = [f1[0]/2.-df1[0]*D+a3*f3[0]/2.-a3*df3[0]*D, f2[0]/2.-df2[0]*D+a4*f4[0]/2.-a4*df4[0]*D, 0, 0]
		# self.C[1,:] = [-f1[0]/4.-a3*f3[0]/4.-D*df1[0]/2.-D*a3*df3[0]/2.-D*df1[1]-D*a3*df3[1], -f2[0]/4.-a4*f4[0]/4.-D*df2[0]/2.-D*a4*df4[0]/2.-D*df2[1]-D*a4*df4[1], 0, 0]
		# self.C[2,:] = [0, 0, f1[1]/2.+df1[1]*D+a3*f3[1]/2.+a3*df3[1]*D, f2[1]/2.+df2[1]*D+a4*f4[1]/2.+a4*df4[1]*D]
		# self.C[3,:] = [0, 0, D*df1[0]+D*a3*df3[0]+f1[1]/4+a3*f3[1]/4.-D*df1[1]/2.-D*a3*df3[1]/2., D*df2[0]+D*a4*df4[0]+f2[1]/4.+a4*f4[1]/4.-D*df2[1]/2.-D*a4*df4[1]/2.]

		# self.S[0] = -flux[n]/2.
		# self.S[1] = (XSf-XSabs)*opt.delta*flux[n]+flux[n]/4.
		# n = opt.numBins-1
		# self.S[2] = -flux[n]/2
		# self.S[3] = (XSf-XSabs)*opt.delta*flux[n]-flux[n]/4.

		# Zero flux boundary condition
		self.C[0,:] = [f1[0]+a3*f3[0], f2[0]+a4*f4[0], 0, 0]
		self.C[1,:] = [df1[0]+a3*df3[0]-df1[1]-a3*df3[1], df2[0]+a4*df4[0]-df2[1]-a4*df4[1], 0, 0]
		self.C[2,:] = [0, 0, f1[1]+a3*f3[1], f2[1]+a4*f4[1]]
		self.C[3,:] = [0, 0, df1[0]+a3*df3[0]-df1[1]-a3*df3[1], df2[0]+a4*df4[0]-df2[1]-a4*df4[1]]
		
		self.S[0] = 0. 
		self.S[1] = (XSf/k-XSr)*opt.delta**2*flux[n]/D

		n = opt.numBins-1
		self.S[2] = 0.
		self.S[3] = (XSf/k-XSr)*opt.delta**2*flux[n]/D

		
###############################################################
#end class 
