# Class containing nuclide objects in a dictionary. 

class Material:
    # Initialization (constructor) routine

    def __init__(self):
        self.name = 'name'

    # Routine for reading material data from a file. 

    def read(self):
        
        
        import numpy as np
        inpFile = open('Materials/MaterialData.inp', 'r')
        self.data={} 
        # Reset n which indicates which line of "Gscat" to read
        n=0
        # Create list of nuclides
        self.materialList = []
        
        
        for line in inpFile:  
            
            
            # Remove trailing whitespace
            line = line.strip()
            # Remove newline characters
            line = line.strip('\n')
            # Remove string after comment character (#)
            line, scratch1, scratch2 = line.partition('#')
            # Skip empty lines  
            if len(line) == 0:
                continue   
                
            
            keyword, arguments = line.split(' ', 1)
            if keyword == 'NumMaterials':
                nMat = int(arguments)

            elif keyword == 'name':
                self.name = arguments
                self.materialList.append(self.name)
                
                
            elif keyword == 'mat':
                self.n = int(arguments)
                
                
            elif keyword == 'numGroups':
                self.nGroups = int(arguments)
                self.Gscat=np.zeros((self.nGroups,self.nGroups))


            elif keyword == 'diffcoef':
                self.D=[]
                for i in range(1,self.nGroups+1):
                    self.D.append(float(line.split(' ',self.nGroups)[i]))
                    
                    
            elif keyword == 'absXS':
                self.absXS=[]
                for i in range(1,self.nGroups+1):
                    self.absXS.append(float(line.split(' ',self.nGroups)[i]))
                    
                    
            elif keyword == 'fisXS':
                self.fisXS=[]
                for i in range(1,self.nGroups+1):
                    self.fisXS.append(float(line.split(' ',self.nGroups)[i]))
                    
                
            elif keyword[:-1] == 'Gscat':  
                for i in range(0, self.nGroups):
                    self.Gscat[n,i]=float(line.split(' ',self.nGroups)[i+1])
                # Next line will read the following line of Gscat
                n=n+1
    
                
            if keyword == 'end':
                n=0
                self.d = {
                    self.name: {
                        'D' : {},
                        'absXS' : {},
                        'fisXS' : {},
                    },
                }
                for i in range(1, self.nGroups+1):
                    self.d[self.name]['D'][i]=self.D[i-1]
                    self.d[self.name]['absXS'][i]=self.absXS[i-1]
                    self.d[self.name]['fisXS'][i]=self.fisXS[i-1]
                    self.d[self.name]['Ex'+str(i)]={}
                    for j in range(1,self.nGroups+1):
                        self.d[self.name]['Ex'+str(i)][j]=self.Gscat[i-1,j-1]
                self.data.update(self.d)
                
                
            else:
                continue
                
                
        # print self.data
        # print self.materialList
