import numpy as np
import math
import matplotlib.pyplot as plt

class DistanceMatrix():
	def __init__(self):
		self.distanceMatrix = []
		self.points = []
		self.periodicMax = None, None, None
		self.distanceDistribution = None
		self.maxShellIndex = 0
		self.distances = None
		self.N = 0

	def setPeriodicBound(self, boxdim):
		# print "Setting periodic max"
		self.periodicMax = boxdim

	def addPoint(self, point):
		# self.points.append(point)
		self.points.append(np.multiply(np.random.rand(3),self.periodicMax))
		self.N = self.N+1
		# print "Adding point: " + str(point)

	def modular1Ddistance(self, x1, x2, boxdim):
		return min([abs(x1-x2), abs(x1-x2+boxdim), abs(x1-x2-boxdim)])

	def distance(self, point1, point2, twoD=False):
		x1, y1, z1 = point1
		x2, y2, z2 = point2
		maxx, maxy, maxz = self.periodicMax
		delx, dely, delz = self.modular1Ddistance(x1,x2, maxx), self.modular1Ddistance(y1,y2, maxy), self.modular1Ddistance(z1,z2,maxz)
		# return (delx**2+dely**2+delz**2)**(.5)
		# only using 2D distance
		# print "Distance: " + str((delx**2+dely**2)**(.5))
		if ((delx**2+dely**2)**(.5))<1:
			if (delx**2+dely**2)**(.5) != 0.0:
				raise Exception("Small Distance: " + str((delx**2+dely**2)**(.5)) + " calculated from: " + str(point1) + " and " + str(point2) + " boxdims is: " + str(self.periodicMax))

		if twoD:
			return (delx**2+dely**2)**(.5)
		else:
			return (delx**2+dely**2+delz**2)**(.5)


	def createDistanceMatrix(self):
		# print "Creating Distance Matrix"
		# print "There are " + str(len(self.points)) + " points."
		# print "Truncating to first 100 points..."
		# self.points = self.points[:100]
		for index1, point1 in enumerate(self.points):
			self.distanceMatrix.append([])
			for index2, point2 in enumerate(self.points):
				
				# if index1<=index2:
				self.distanceMatrix[index1].append(self.distance(point1,point2, twoD=True))
				# Using D_ij = D_ji 
				# if index1>index2:
				# 	self.distanceMatrix[index1].append(self.distanceMatrix[index2][index1])
		self.distanceMatrix = np.asarray(self.distanceMatrix)

	def calculateDistanceDistribution(self):
		# print "Calculating Distance Distribution"
		
		self.createDistanceMatrix()
		# print max(self.distanceMatrix)
		distances = []
		# print "Total rows: " + str(len(self.distanceMatrix))
		for row in self.distanceMatrix:
			# print "Row number: " + str(row)
			for distance in row:
				if distance == 0:
					continue
				distances.append(distance)
		self.distances = distances
		distances.sort()
		print distances[:3]

	def pyplotDistanceDistribution(self, numbins=100):
		n, bins, patches = plt.hist(self.distances, bins = numbins)
		plt.clf()
		dr = bins[1]-bins[0]
		x, y, z = self.periodicMax
		density = self.N/(x*y*z)
		# weights = [1./(2*math.pi*r*dr*self.periodicMax[2]*self.N*density) for r in self.distances]
		# weights = np.ones(bins)
		# print self.distances
		return np.histogram(self.distances, bins=numbins)

	def showDistanceDistribution(self, numbins=100):
		n, bins, patches = self.pyplotDistanceDistribution(numbins = numbins)
		plt.plot(bins[1:], n)
		plt.plot(bins[1:], [1.]*(len(bins)-1), "--", color = "gray")
		# print n
		print "box is : " + str(self.periodicMax)
		plt.ylim(0,2)
		plt.xlim(0, min(self.periodicMax[0],self.periodicMax[1])/2)
		plt.title("g(r)")
		plt.show()
		plt.clf()

# dm = DistanceMatrix()
# N = 500
# boxx, boxy, boxz = 10., 10., 10.
# dm.setPeriodicBound([boxx,boxy,boxz])
# for n in range(N):
# 	dm.addPoint(np.multiply(np.random.rand(3),dm.periodicMax))
# dm.calculateDistanceDistribution()
# dm.showDistanceDistribution()

# dm.addPoint((1,2,3))
# dm.addPoint((9,2,3))
# dm.addPoint((5,3,1))
# dm.createDistanceMatrix()
# # print str(dm.distanceMatrix)
# dm.calculateDistanceDistribution()

