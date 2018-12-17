# Reads input
# I chose to use number of bins rather than the mesh size to ensure that a user couldn't use 
# a mesh size that didn't fit into the length. But the capability to use meshSize instead of 
# numBins is here. 

class DiffusionOpts1D:
    
    def __init__(self):
        self.length=-1
        
    def read(self, filename):
        inpFile = open(filename, 'r')
        
        
        
        for line in inpFile:  


            #Remove trailing white space
            line = line.strip()
            #Remove newline characters
            line = line.strip('\n')
            #Remove string after comment character (#)
            line, scratch1, scratch2 = line.partition('#')
            #Skipe empty lines
            if len(line) == 0:
                continue


            keyword, arguments = line.split(' ', 1)
            if keyword == 'length':
                self.length = float(arguments)
                
            elif keyword == 'numgroups':
                self.numGroups = int(arguments)
                
            elif keyword == 'numBins':
                self.numBins = int(arguments)

            elif keyword == 'meshSize':
                self.delta = float(arguments)

            elif keyword == 'FisConvergeError':
                self.FisConvError = float(arguments)

            elif keyword == 'FluxConvergeError':
                self.FluxConvError = float(arguments)

            elif keyword == 'problem':
                self.pb_num = int(arguments)
                
            else:
                continue

        self.delta = self.length/self.numBins
        # self.numBins = int(self.length/self.delta)