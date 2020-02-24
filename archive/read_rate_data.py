# Class for reading and storing rate data from different locations
# Created based on read_EIRENE.py by holm10 on Jan 27 2020
# Changelog:
#   200127 - Rewrote read_EIRENE into a class



class RATE_DATA:
    def __init__(self,amjuel='amjuel.tex',hydhel='hydhel.tex',h2vibr='h2vibr.tex',ADAS='ic#h0-1.dat',UE='ehr1.dat',path='.'):
        self.reactions= {   'AMJUEL': {'settings' : [500, ['b0','0','a0','h0','p0','k0'], ['a0','h0','p0','k0'], 45], 'path' : amjuel },
                            'HYDHEL': {'settings' : [150, ['b0','0','a0','h0'], ['a0','h0'], 80], 'path' : hydhel },
                            'H2VIBR': {'settings' : [0, ['b0','0','a0','h0'], ['a0','h0'], 20], 'path' : h2vibr },
                            'UE': {},
                            'ADAS': {},
                        }

        for rate in ['AMJUEL','H2VIBR','HYDHEL']:
            self.read_EIRENE(self.reactions[rate]['path'],self.reactions[rate],self.reactions[rate]['settings'],path=path)
        self.read_ADAS(ADAS,self.reactions['ADAS'])
        self.read_UE(UE,self.reactions['UE'])




    def read_EIRENE(self,fname,reactions,settings,path='.'):
        from numpy import zeros
        from csv import reader
        # Script parsing Hydhel.tex
        lst,book=[],[]
        # Open reader object
        rd=reader(open('{}/{}'.format(path,fname),'rt'),delimiter=' ')
        # Store book in list
        for row in rd:
            lst.append(row)
        # Strip all empty arrays
        lst=list(filter(None,lst))

        for row in lst:
            book.append(list(filter(None,row)))

        i=settings[0] # Read from top
        while True:
            # Loop through rows looking for reactions
            if book[i][0]=="Reaction":
                reaction=book[i][1] # Store reaction name
                i+=1
                # Loop through reaction looking for coeffs
                while True:
                    if book[i][0] not in settings[1]:
                        i+=1
                     #   break
                    # Break if wrong fits
                    elif book[i][0] in settings[2]:
                        break
                    # We are in a T-fit
                    elif book[i][0]=='b0':
                        coeff=[]
                        # Parse the next three lines
                        for j in range(3):
                            for k in range(3):
                                coeff.append(float(book[i+j][1+k*2].replace(',','').replace('D','E')))
                        i+=1
                        reactions[reaction]=coeff
                        break
                    # Wea re in a (T,E)-fit
                    elif book[i][0]=='0':
                        coeff=zeros((9,9))
                        # Parse the whole data block in one
                        for j in range(3):
                            for k in range(9):
                                for l in range(3):
                                    coeff[k,j*3+l]=float(book[i+k+j*9+j*2][l+1].replace('D','E'))
                        # Store the coefficients
                        # TODO: figure out better way to kill off ne,T fits??
                        if reaction not in ['2.2.14','2.0l2']:
                            reactions[reaction]=coeff
                        
                        i+=9+3*3+2*2 # Set counter after block
                        break
                    if i>=len(book)-settings[3]: # Omit last lines 
                        break

            i+=1
            if i>=len(book)-settings[3]: # Omit last lines
                break


    def read_ADAS(self,fname,reactions,path='.'):
        # Constants turning the values into 
        cm1=1/8065.6
        kB=8.621738e-5
        # The ADAS reading routine follows the ADAS documentation: only parameters relevant to the CRM are extracted
        with open('{}/{}'.format(path,fname)) as f:  # Open the file
            l=f.readline() # Discard data
            l=f.readline().split() # Read first data line
            # Loop through energy level data
            while l[0]!='-1':
                reactions['E'+l[0]]=float(l[4])*cm1 # Store as eV in reactions['E']
                l=f.readline().split()
            # Read and store temperature point data
            l=f.readline().split()
            reactions['T']=[float(x.replace('+','e+').replace('-','e-'))*kB for x in l[2:]]
            # Read first rate line
            l=f.readline()
            # Read all unspecified data
            while l[0]==' ':
                if l.strip()=='-1': break
                [ul,ll]=l[:8].split()
                l=l[8:]
                l=[l[i:i+8].strip() for i in range(0,len(l),8)][:-1]
                # Store excitation and relaxation data
                reactions[ul+'-'+ll]=float(l[0][:-3]+'e'+l[0][-3:])
                reactions[ll+'-'+ul]=[float(x[:-3]+'e'+x[-3:]) for x in l[1:len(reactions['T'])+1]]
                l=f.readline()


    def read_UE(self,fname,reactions,path='.'):
        from numpy import array,transpose
        with open('{}/{}'.format(path,fname)) as f:  # Open the file
            # Read the ionizationa nd recombination data (1st and 2nd blocks in ehr1.dat)
            for k in ['IONIZ','REC']:
                reactions[k]=[]
                l=f.readline() # Discard header
                # Loop through data blocks
                for i in range(15): # Density points
                    buff=[]    # Buffer
                    for j in range(12): # Temperature points
                        l=f.readline().split() 
                        # Arrange temperature points into 1D array
                        if j not in [0,11]:
                            for x in l:
                                buff.append(float(x))
                    # Append temperature points to reactions['IONIZ/REC'] for each density point
                    reactions[k].append(buff)
                reactions[k]=transpose(array(reactions[k])) # Convert 2D list into array(T,n)

    def get_coeff(self,database,reaction):
        return(self.reactions[database][reaction])



                



