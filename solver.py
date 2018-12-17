# class to solve by power iteration

import numpy as np

class Solve:
	
	def __init__(self, opt):
		self.J = np.zeros(opt.numBins+1)
		self.flux = np.zeros(opt.numBins*opt.numGroups)

###############################################################
		
	def solve(self, opt, Ainv, F):
		# Function to solve for flux

		# Create local array variables.
		
		RMSflux = np.zeros(opt.numBins*opt.numGroups)
		self.B = np.zeros(opt.numBins*opt.numGroups)
		RMSsource = np.zeros(opt.numBins*opt.numGroups)

		# Initialize local variables
		ERRflux = 1000
		ERRsource = 1000
		j = 0
		# Initial eigenvalue guess
		self.k = 1
		# Initial flux guess
		self.flux[:] = 1/opt.length
		self.B[:] = 1

		# while j < 10:

		while ERRflux > opt.FluxConvError or ERRsource > opt.FisConvError:

			# Reset local variables
			j = j+1
			RMSsource[:] = 0
			n = 0
			B = []

			# Save previous source so that we can compute the
			# residual error for the iteration
			lastFlux = self.flux
			lastB = self.B

			self.B = np.dot(F,self.flux)/self.k
			
			self.flux = np.dot(Ainv,self.B)

			# Calculate the fission source in each spatial bin
			self.k = sum(np.dot(F,self.flux))/sum(self.B) #sum(np.dot(F,lastFlux))*k

			# Perform the source normalization
			self.flux[:] = self.flux[:]/np.linalg.norm(self.flux)

			# Calculate the relative difference in the source between
			# consecutive iterations and take the infinity norm.
			RMSflux[:] = abs((lastFlux[:]-self.flux[:])/self.flux[:])
			for i in range(0,len(self.B)):
				if self.B[i] != 0:
					n = n+1
					B.append(self.B[i])
					RMSsource[i] = abs((lastB[i]-self.B[i])/self.B[i])

			ERRflux = max(RMSflux)/len(self.flux)
			ERRsource = max(RMSsource)/n
			# print(j)
			# print("Error:", ERRflux, ERRsource)

			# Print statement to show eigenvalue convergence by iteration
			# print("Eigenvalue: ", self.k, "(", ERRflux, ";", ERRsource, ") Iteration:", j)
			# print(self.flux)
		# print("Eigenvalue: ", self.k, "(", ERRflux, ";", ERRsource, ") Iteration:", j)
		self.it = j
		# print j
		# avgB = np.average(B)
		# self.peak = max(self.B)/np.average(B)
		# self.peakLoc = abs(opt.length/2.-B.index(max(B))*opt.delta)
		

###############################################################

	def coefs(self, Cinv, S):
		# Function to calculate flux expansion coefs
		self.a = np.dot(Cinv,S)
		# print(self.a)

###############################################################

	def current(self, n, a, opt, D, XSr, flux, Matrix, 
		E, f1, df1, f2, df2, f3, df3, f4, df4):
	    # Function to calculate current from flux expansion coefs

		XSd = D/opt.delta**2

		a1 = a[2]
		a2 = a[3]
		a3 = -a1*5.*XSr/(3.*(60.*XSd+XSr))
		a4 = a2*7*XSr/(420.*(XSd+3.*XSr))

		if n == 0:

			# Zero flux boundary condition

			a1 = a[0]
			a2 = a[1]
			a3 = -a1*5.*XSr/(3.*(60.*XSd+XSr))
			a4 = a2*7*XSr/(420.*(XSd+3.*XSr))

			self.J[n] = -D/opt.delta*(a1*df1[0]+a2*df2[0]+a3*df3[0]+a4*df4[0])
			# print(a1*df1[0]+a2*df2[0]+a3*df3[0]+a4*df4[0])
			# print(a1*df1[1]+a2*df2[1]+a3*df3[1]+a4*df4[1])

			a1 = a[2]
			a2 = a[3]
			a3 = -a1*5.*XSr/(3.*(60.*XSd+XSr))
			a4 = a2*7*XSr/(420.*(XSd+3.*XSr))
			# print(a1*df1[1]+a2*df2[1]+a3*df3[1]+a4*df4[1])
			# print(a1*df1[0]+a2*df2[0]+a3*df3[0]+a4*df4[0])

			self.J[opt.numBins] = -D/opt.delta*(a1*df1[1]+a2*df2[1]+a3*df3[1]+a4*df4[1])

		else:
			self.J[n] = -D/opt.delta*(a1*df1[0]+a2*df2[0]+a3*df3[0]+a4*df4[0])
			Matrix.Dh[n] = (-Matrix.Ds[n]*(flux[n]-flux[n-1])-self.J[n])/(flux[n]+flux[n-1])

###############################################################

	def node_bal_check(self, delta, flux, XSr, XSf, A, Ds, Dh, k):
		# Function to check nodal balance in each cell
		self.ERR = []
		self.BoundERR = []
		self.AERR = []
		self.CMFDerr = []

		for i in range(0,len(self.J)-1):
			if i >= 1 and i < len(self.J)-2:
				Jr = -Ds[i+1]*(flux[i+1]-flux[i])-Dh[i+1]*(flux[i+1]+flux[i])
				Jl = -Ds[i]*(flux[i]-flux[i-1])-Dh[i]*(flux[i-1]+flux[i])
				self.CMFDerr.append(abs(Jr-Jl+(XSr-XSf/k)*delta*flux[i]))
				# print Jl
				if i >= 2:
					self.AERR.append(abs(A[i,0]-A[i-1,2]))
				self.ERR.append(abs(self.J[i+1]-self.J[i]+(XSr-XSf/k)*delta*flux[i]))
		self.BoundERR.append(abs(self.J[i+1]-self.J[i]+(XSr-XSf/k)*delta*flux[i]))

		# print("CMFD balance is verified to be %f" %(max(self.CMFDerr)))
		print("NEM balance is verified to be %f and coef error is %f" %(
			max(self.ERR), max(self.AERR)))
		print("But at the boundaries, nodal balance is not satisfied: %f" %(max(self.BoundERR)))

		# print self.J[20]-self.J[19]+(XSr-XSf/k)*delta*flux[19]

		# print counter
		# print(self.ERR)


###############################################################
#end class
		








