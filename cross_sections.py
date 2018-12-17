# Class to fill cross section vectors

import numpy as np

class CrossSections:

###############################################################

	def __init__(self, opts):
		self.abs = []
		self.fis = []
		self.removal = []
		self.D = []
		self.Gscat = np.zeros((opts.numBins,opts.numGroups*opts.numGroups))

###############################################################

	def problem1(self, opts, data):

		nBins = opts.numBins
		nGrps = opts.numGroups

		for g in range(1,nGrps+1):
			# If the outgoing group is 1, the inscatter group is 2, and vice versa.
			# This strucutre allows for upscatter.
			if g == 1:
				gIn = 2
			else:
				gIn =1
			for i in range(nBins*(g-1),nBins*g):
				self.fis.append(data['fuel1']['fisXS'][g])
				self.abs.append(data['fuel1']['absXS'][g])
				self.removal.append(data['fuel1']['absXS'][g]+data['fuel1']['Ex'+str(g)][gIn])
				self.D.append(data['fuel1']['D'][g])


		# Gscat matrix
		for i in range(0,nBins):
			count = 1
			gcount = 1
			for g in range(0,nGrps*nGrps):
				self.Gscat[i,g] = data['fuel1']['Ex'+str(gcount)][count]
				count = count+1
				if count == nGrps+1:
					count = 1
					gcount = gcount+1

###############################################################

	def problem2(self, opts, data):

		nBins = opts.numBins
		nGrps = opts.numGroups

		for g in range(1,nGrps+1):
			if g == 1:
				gIn = 2
			else:
				gIn = 1
			for i in range(nBins*(g-1),nBins*g):
				if (i < nBins*(g-1)+5) or (i >=nBins*g-5):
					self.fis.append(data['baff_reflect']['fisXS'][g])
					self.abs.append(data['baff_reflect']['absXS'][g])
					self.removal.append(data['baff_reflect']['absXS'][g]+data['baff_reflect']['Ex'+str(g)][gIn])
					self.D.append(data['baff_reflect']['D'][g])
				else:
					self.fis.append(data['fuel1']['fisXS'][g])
					self.abs.append(data['fuel1']['absXS'][g])
					self.removal.append(data['fuel1']['absXS'][g]+data['fuel1']['Ex'+str(g)][gIn])
					self.D.append(data['fuel1']['D'][g])


		# Gscat matrix
		for i in range(0,nBins):
			count = 1
			gcount = 1
			for g in range(0,nGrps*nGrps):
				if i < 5 or i >= nBins-5:
					self.Gscat[i,g] = data['baff_reflect']['Ex'+str(gcount)][count]
				else:
					self.Gscat[i,g] = data['fuel1']['Ex'+str(gcount)][count]
				count = count+1
				if count == nGrps+1:
					count = 1
					gcount = gcount+1

###############################################################

	def problem3(self, opts, data):

		nBins = opts.numBins
		nGrps = opts.numGroups

		for g in range(1,nGrps+1):
			if g == 1:
				gIn = 2
			else:
				gIn = 1
			for i in range(nBins*(g-1),nBins*g):
				if (i < nBins*(g-1)+25) or (i >=nBins*g-25):
					self.fis.append(data['baff_reflect']['fisXS'][g])
					self.abs.append(data['baff_reflect']['absXS'][g])
					self.removal.append(data['baff_reflect']['absXS'][g]+data['baff_reflect']['Ex'+str(g)][gIn])
					self.D.append(data['baff_reflect']['D'][g])
				else:
					self.fis.append(data['fuel1']['fisXS'][g])
					self.abs.append(data['fuel1']['absXS'][g])
					self.removal.append(data['fuel1']['absXS'][g]+data['fuel1']['Ex'+str(g)][gIn])
					self.D.append(data['fuel1']['D'][g])


		# Gscat matrix
		for i in range(0,nBins):
			count = 1
			gcount = 1
			for g in range(0,nGrps*nGrps):
				if i < 25 or i >= nBins-25:
					self.Gscat[i,g] = data['baff_reflect']['Ex'+str(gcount)][count]
				else:
					self.Gscat[i,g] = data['fuel1']['Ex'+str(gcount)][count]
				count = count+1
				if count == nGrps+1:
					count = 1
					gcount = gcount+1

###############################################################

	def problem4(self, opts, data):

		nBins = opts.numBins
		nGrps = opts.numGroups

		for g in range(1,nGrps+1):
			if g == 1:
				gIn = 2
			else:
				gIn = 1
			for i in range(nBins*(g-1),nBins*g):
				if (i < nBins*(g-1)+25) or (i >=nBins*g-25):
					self.fis.append(data['baff_reflect']['fisXS'][g])
					self.abs.append(data['baff_reflect']['absXS'][g])
					self.removal.append(data['baff_reflect']['absXS'][g]+data['baff_reflect']['Ex'+str(g)][gIn])
					self.D.append(data['baff_reflect']['D'][g])
				elif (i < nBins*(g-1)+40) or (i >=nBins*g-40):
					self.fis.append(data['fuel4']['fisXS'][g])
					self.abs.append(data['fuel4']['absXS'][g])
					self.removal.append(data['fuel4']['absXS'][g]+data['fuel4']['Ex'+str(g)][gIn])
					self.D.append(data['fuel4']['D'][g])
				else:
					self.fis.append(data['fuel3']['fisXS'][g])
					self.abs.append(data['fuel3']['absXS'][g])
					self.removal.append(data['fuel3']['absXS'][g]+data['fuel3']['Ex'+str(g)][gIn])
					self.D.append(data['fuel3']['D'][g])


		# Gscat matrix
		for i in range(0,nBins):
			count = 1
			gcount = 1
			for g in range(0,nGrps*nGrps):
				if i < 25 or i >= nBins-25:
					self.Gscat[i,g] = data['baff_reflect']['Ex'+str(gcount)][count]
				elif i < 40 or i >= nBins-40:
					self.Gscat[i,g] = data['fuel4']['Ex'+str(gcount)][count]
				else:
					self.Gscat[i,g] = data['fuel3']['Ex'+str(gcount)][count]
				count = count+1
				if count == nGrps+1:
					count = 1
					gcount = gcount+1


###############################################################

	def problem5(self, opts, data):

		nBins = opts.numBins
		nGrps = opts.numGroups

		for g in range(1,nGrps+1):
			if g == 1:
				gIn = 2
			else:
				gIn = 1
			for i in range(nBins*(g-1),nBins*g):
				if (i < nBins*(g-1)+23) or (i >=nBins*g-23):
					self.fis.append(data['reflector']['fisXS'][g])
					self.abs.append(data['reflector']['absXS'][g])
					self.removal.append(data['reflector']['absXS'][g]+data['reflector']['Ex'+str(g)][gIn])
					self.D.append(data['reflector']['D'][g])
				elif (i < nBins*(g-1)+25) or (i >=nBins*g-25):
					self.fis.append(data['baffle']['fisXS'][g])
					self.abs.append(data['baffle']['absXS'][g])
					self.removal.append(data['baffle']['absXS'][g]+data['baffle']['Ex'+str(g)][gIn])
					self.D.append(data['baffle']['D'][g])
				elif (i < nBins*(g-1)+40) or (i >=nBins*g-40):
					self.fis.append(data['fuel4']['fisXS'][g])
					self.abs.append(data['fuel4']['absXS'][g])
					self.removal.append(data['fuel4']['absXS'][g]+data['fuel4']['Ex'+str(g)][gIn])
					self.D.append(data['fuel4']['D'][g])
				else:
					self.fis.append(data['fuel3']['fisXS'][g])
					self.abs.append(data['fuel3']['absXS'][g])
					self.removal.append(data['fuel3']['absXS'][g]+data['fuel3']['Ex'+str(g)][gIn])
					self.D.append(data['fuel3']['D'][g])

		# Gscat matrix
		for i in range(0,nBins):
			count = 1
			gcount = 1
			for g in range(0,nGrps*nGrps):
				if i < 23 or i >= nBins-23:
					self.Gscat[i,g] = data['reflector']['Ex'+str(gcount)][count]
				elif i < 25 or i >= nBins-25:
					self.Gscat[i,g] = data['baffle']['Ex'+str(gcount)][count]
				elif i < 40 or i >= nBins-40:
					self.Gscat[i,g] = data['fuel4']['Ex'+str(gcount)][count]
				else:
					self.Gscat[i,g] = data['fuel3']['Ex'+str(gcount)][count]
				count = count+1
				if count == nGrps+1:
					count = 1
					gcount = gcount+1


###############################################################


#end class


