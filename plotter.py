#class to plot


class Plotter:

####################################################################
    def __init__(self):
    	self


####################################################################

    def plot_flux(self, opt, flux, j):

    	import numpy as np
    	import matplotlib.pyplot as plt

        plt.clf()
        nBins = opt.numBins
        nGrps = opt.numGroups

        s = np.zeros(nBins+1)
        m = np.zeros(nBins)

        for n in range(0,nBins+1):
            s[n] = (n)*opt.delta
        for n in range(0,nBins):
            m[n] = s[n] + opt.delta/2.

    	groups = []
    	for g in range(1,nGrps+1):
    		l = 0
    		x_group = np.zeros(nBins)
    		for i in range(nBins*(g-1),nBins*g):
    			x_group[l] = flux[i]
    			l = l+1


    		groups.append('group%i' %g)
    		colors = ['r','b','g','o','n','p']
    		plt.plot(m, x_group, colors[g-1])
    		plt.ylabel('flux')
    		plt.savefig('./flux')
            # plt.savefig('./flux%i' %j)

####################################################################

    def plot_flux_change(self, opt, First, Last):

        import numpy as np
        import matplotlib.pyplot as plt

        plt.clf()
        nBins = opt.numBins

        s = np.zeros(nBins+1)
        m = np.zeros(nBins)

        for n in range(0,nBins+1):
            s[n] = (n)*opt.delta
        for n in range(0,nBins):
            m[n] = s[n] + opt.delta/2.

        plt.plot(m, abs(First[:]-Last[:])/Last[:])
        plt.ylabel('Relative change in flux from first iteration')
        plt.savefig('./fluxChange')

####################################################################

    def plot_node_err(self, ERR, BoundERR, AERR):

        import numpy as np
        import matplotlib.pyplot as plt

        plt.clf()
        plt.plot(ERR)
        plt.title('Node balance error in cells')
        plt.savefig('./error1')
        plt.clf()
        plt.plot(BoundERR)
        plt.title('Node balance error in boundary cells')
        plt.savefig('./error2')
        plt.clf()
        plt.plot(AERR)
        plt.title('Coefficient error')
        plt.savefig('./error3')
        plt.clf()


####################################################################
#end class


