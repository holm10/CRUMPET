# CRUM collisional-radiative model class reaction.py
# Separated from CRUM.py by holm10
# Changelog
# 200205 - Separated from CRUM.py #holm10
from CRUM.tools import TOOLS


class CRM(TOOLS):


    def __init__(self,species,bg,reactionlist,verbose,NP,path='.',
                                            vmax=14,nmax=8,rdata=None):
        ''' Creates a CRM class, at the heart of CRUM
            __init__(species,reactions,settings)

            species     -   List of strings with each species handle/identifier
            reactions   -   List of reaction objects to be included in the CRM
            settings    -   List of paramters controlling CRM behaviour
                            settings[0] - Verbose, switch to show more output
                            settings[1] - Np, number of P-space species taken as the first 
                                          members of species
                            settings[2] - List of same length as species defining the 
                                          initial state of the CRM species

            Optional parameters
            path ('.')  -   Path to CRUm run directory
        '''
        from os import mkdir,getcwd
        from datetime import datetime

        # Store class objects
        self.species=species
        self.bg=bg
        self.slist=list(self.species)
        self.verbose=verbose
        self.Np=NP
        self.path=path
        self.reactions=[]
        self.vmax=vmax
        self.nmax=nmax


        self.setup_reactions(reactionlist,rdata)

        # Define reactions for UEDGE raditation
        #self.ionizrad=REACTION('IONIZRAD','UE',
                #self.get_coeff('UE','IONIZRAD'),'UE',['','',''],bg,species,
                #['','',None,None,[0,0,0,0]])
        #self.recrad=REACTION('RECRAD','UE',self.get_coeff('UE','RECRAD'),
                #'UE',['','',''],bg,species,['','',None,None,[0,0,0,0]])



        # Ensure that there is a logs directory under the run path
        try:
            mkdir('{}/logs'.format(self.path))
        except:
            pass


        # Write a log of the CRM setup path/logs
        with open('{}/logs/setup.log'.format(self.path),'w') as f:
            f.write('CRUM run in {} on {}\n'.format(getcwd(),
                                            str(datetime.now())[:-7]))
            f.write('Defined species:\n')
            for i in self.slist:
                f.write('    {}\n'.format(i))

            f.write('Defined reactions:\n')
            for r in self.reactions:
                f.write('{}\n'.format(r.print_reaction()))
        # Output to stdout if run verbosely
        if self.verbose:
            with open('logs/setup.log','rt') as f:
                for l in f:
                    print(l.strip())
        # Do the same for a Diagnostic rate matrix displaying 
        # reaction correlations
        self.DIAGNOSTIC()



# %%%%%%%%%%%%%% INPUT FILE READING TOOLS %%%%%%%%%%%%%%%%%

    def setup_ADAS(self,x,ID,rdata,rlist,plist,K):
        from CRUM.reactions import REACTION
        if ID=='EXCITATION': # Electron impact excitation
            rn=range(x+1,self.nmax+1)   # Excitation only possible between
                                        # current state and nmax
            fit='ADAS'          # ADAS-type fit
            Tarr=rdata[database]['T']  # Temperature array for interpolation


        elif ID=='RELAXATION': # Radiative relaxation
            rn=range(1,x)   # Relaxation only possible to lower states
            fit='COEFFICIENT' # Coefficient-like decay
            Tarr=None   # Not T-dependent

        # Loop through each of the available final states
        for y in rn:
    
            try:
                _name='{}_{}-{}'.format(ID,x,y)
                _ratecoeff=rdata['ADAS']['{}-{}'.format(x,y)] #Get coefficients
                self.reactions.append(REACTION(_name,'ADAS',_ratecoeff,fit,
                                            [rlist,plist,K],bg,species,Tarr))

                
            except:
                pass


    def setup_APID(self,x,ID,rlist,plist,K):
        from CRUM.reactions import REACTION

        self.reactions.append(REACTION(self.XY2num(ID,x),'APID',
                        x,'APID',[rlist,plist,K],self.bg,self.species))

    def setup_Johnson(self,i,ID,r):
        ''' Johnson's approximation of Einstein coefficients '''
        from CRUM.reactions import REACTION
        from numpy import pi

        def g(i,f):
            g=[ 1.133*(f==1) + 1.0785*(f==2) + (0.9935 + 0.2328/f \
                                                    - 0.1296/f**2)*(f>2),
                    -0.4059*(f==1) -0.2319*(f==2) - ((0.6282 - 0.5598/f \
                                                    + 0.5299/f**2)/f)*(f>2),
                    0.07014*(f==1) + 0.02947*(f==2) + ((0.3887 - 1.181/f \
                                                    + 1.470/f**2)/f**2)*(f>2) ]
            x=1-(f/i)**2
            return g[0] + g[1]/x + g[2]/x**2

        h=1.054571596e-27
        c=2.99792458e10
        me=9.10938188e-28
        e=4.80274e-10
        I=(me * e**4) / (2 * h**2)

        for f in range(1,i):
            res= (2**6 *e**2 * I**2) / (3**1.5 * pi * me *c**3 * h**2)
            freq=(1/f**2 - 1/i**2)
            Afac=(res*g(i,f))/(freq*(i**5)*(f**3))

            rlist=[self.XY2num(a,i,f) for a in r[1]]  # Reactants 
            plist=[self.XY2num(a,i,f) for a in r[2]]  # Fragments


            _rlist=[rlist,plist,r[-1]]
            self.reactions.append(REACTION(self.XY2num(ID,i,f),'JOHNSON',Afac,
                                    'COEFFICIENT',_rlist,self.bg,self.species))






    def setup_reactions(self,reactionlist,rdata):
        from CRUM.reactions import REACTION
        from numpy import zeros,pi
        
        ''' LOOP OVER ALL THE DEFINED REACTIONS '''
        for r in reactionlist:
            # Split the definitions into names and databases
            database=r[0].split('_')[0].upper()
            try:
                ID=r[0].split('_')[1].upper()
                spec=''.join(r[1]+r[2]).upper()
            except:
                ID='CUSTOM'
                spec=''

            # Loop through states, if necessary. Dynamicallt set boundaries
            # according to electronic or vibrational transitions
            for x in range('N=$' in spec,
                        1+('V=$' in spec)*self.vmax+('N=$' in spec)*self.nmax):

                # Vibrational/electronic dependence present
                if '$' in ID:

                    # Substitute state into reactants and product strings 
                    rlist=[self.XY2num(i,x) for i in r[1]]  # Reactants 
                    plist=[self.XY2num(i,x) for i in r[2]]  # Fragments

                    # %%% ADAS rates detected %%%
                    if database=='ADAS':
                        self.setup_ADAS(x,ID,rdata,rlist,plist,r[-1])
                      
                    # %%% APID rates detected %%%
                    elif database=='APID':
                        self.setup_APID(x,ID,rlist,plist,r[-1])
                        
          
                    # %%% APID rates detected %%%
                    elif database=='JOHNSON':
                        self.setup_Johnson(x,ID,r)
                      
                        
                    # %%% Neither of the above: rate for transitions %%%
                    else:
                        # Assume ladder-like vibrational transitions (+/-1) only
                        for y in range(-1*('&' in ID),2,2):

                            # Limit transitions to [0,vmax]
                            if x+y in range(self.vmax+1):
                                # Retain intial and final states in name
                                vID=self.XY2num(ID,x,x+y) 
                                # List of reactants with initial states
                                rlist=[self.XY2num(i,x,x+y) for \
                                                            i in r[1]]
                                # List of products with final states
                                plist=[self.XY2num(i,x,x+y) for \
                                                            i in r[2]]
                                _name=vID
                                _database=database
                                _ratecoeff=rdata[database][vID] 
                                _rtype='RATE'
                                _rlist=[rlist,plist,r[-1]]

                                self.reactions.append(REACTION(
                                            _name,_database,_ratecoeff,_rtype,
                                            _rlist,self.bg,self.species)      )
                        
                #%%%%% Read custom rates %%%%%
                elif database=='CUSTOM':
                    self.setup_custom(r[1].strip(),database)

                #%%% EIRENE/UEDGE-DEGAS rates %%%
                elif database in ['HYDHEL','AMJUEL','H2VIBR','UE']:
                    self.reactions.append(REACTION(ID,database,
                                rdata[database][ID],
                                'RATE'*(database!='UE')+'UE'*(database=='UE'),
                                [r[1],r[2],r[-1]],self.bg,self.species))

                #%%% Fell through loop %%%
                else:
                    print('''Database "{}" not recognized! 
                                Aborting.'''.format(database))
                    return







    def setup_custom(self,fname,database):
        from CRUM.reactions import REACTION
        from numpy import zeros,pi
        
        # Parse the custom rate file into a list and retain subcards
        data,_,subcards=self.file2list(self.path,fname)
        _database=data[0].strip() # Database is defined at fist line

        # Loop through the rates, each a separate subcard
        for i in range(len(subcards)-1):
            subc=subcards[i] # Helper index
            fit=data[subc].split()[1].upper() # Rate typy, first subcard entry
            name=data[subc].split()[2] # Name/ID of reacrtion - second entry
            reactants,fragments=data[subc+1].split(' > ') # Reactants and 
                                    # fragments are defined on the second line
            K='0' # Kinetic energy exchange in reaction is 0 unless defined

            # Execute if no vib dependence, loop if vibr. dep. process
            for j in range(1+('$' in name)*self.vmax):

                # %%% Vibrationally dependent process %%%
                if '$' in name:

                    # Read kinetic energy for each process
                    for k in range(0,100):
                        if data[k][0]=='K': 
                            K=data[k].strip().split('=')[-1]
                        elif data[k][0]=='v':
                            m=k
                            break

                    # Write data
                    _name=self.XY2num(name,j)
                    rlist=self.XY2num(reactants,j).strip().split(' + ') 
                                                # Reactants w/ v-level number
                    plist=self.XY2num(fragments,j).strip().split(' + ') 
                                                # Products w/ v-level number
                    _ratecoeff=[float(x) for x in data[subc+m+j*2].split()] 
                                                # Coefficients for v-level

                # %%% Specified rate %%%
                else:
                    
                    # Read the kinetic energy of the process
                    for k in range(2,100):
                        if data[subc+k][0]=='K': 
                            K=data[subc+m].strip().split('=')[-1]
                        else:
                            m=k
                            break
        
                    # Cross-section as defined in SAWADA 95 has special form
                    if fit=='SIGMA': 
                        _ratecoeff=[float(x) for x in data[subc+m].split()] 
                    # Other processes have pre-defined form
                    else: _ratecoeff=float(data[subc+m])
                    
                    # Write data
                    _name=name
                    rlist=reactants.strip().split(' + ') # Reactants
                    plist=fragments.strip().split(' + ') # Fragments

                # Store reaction
                self.reactions.append(REACTION(_name,_database,_ratecoeff,
                                    fit,[rlist,plist,K],self.bg,self.species))








            





# %%%%%%%%%%%%%% ENF OF INPUT FILE READING TOOLS %%%%%%%%%%%%%%%%%










        """

    def full_nt(self,Te,ne,t,Ti=None,ni=None,E=0.1,n=None,Sext=True):
        ''' Calculates the full CRM density evolution in a 1cm**3 box up to t
            gl_nt(Te,ne,t)
            Te  -   electron background temperature in box [eV]
            ne  -   electron background density in box [cm**-3]
            t   -   final time of evolution [s]

            Optional parameters
            Ti (None)   -   ion background temperature in box (=Te if None) [eV]
            ni (None)   -   ion background density in box (=ne if None) [cm**-1]
            E (0.1)     -   target particle energy [eV]
            n (None)    -   initial species distribution (=n0 specified in input if None)
                            Array of same length as species, order according to 'SPECIES' card
            Sext (True) -   Include external source (from background plasma reactions into CRM species)
            
            Uses the CRM function full_nt
        '''
        return self.full_nt(Te,ne,t,Ti,ni,E,n,Sext)
    """


    def populate(self,mode,Te,ne,Ti=None,ni=None,E=0,rad=True,Sind=None,Tm=False,Ton=True,Iind=0):
        ''' Function populating a matrix according to the chosen mode 
            populate(mode,Te,ne,*keys)
            mode    -   Matrix writing mode
                            'diagnostic'    -   Creates a 2D list of reactions handles
                            'R'             -   Creates a matrix of rate coefficients (cm**3/s)
                            'M'             -   Creates a matrix of rates (s**-1)
            Te      -   Background plasma electron temperature [eV]
            ne      -   Background plasma electron density [cm**-3]

            Optional parameters
            Ti (None)   -   Background plasma ion temperature [eV]. Ti=Te assumed if None
            ni (None)   -   Background plasma ion density [cm**-3]. ni=ne assumed if None
            E (0.1)     -   Target particle energy [eV]
            rad

            Returns
            matrix,ext_source
    
            matrix      -   The matrix with the requested elements
            ext_source  -   Vector containing the external source contributions
                            to each species. Note that the external sources always are
                            returned as volumetric rates (cm**-3 s**-1), which differ from the
                            rate matrix
            
        '''
        from numpy import zeros,array,sum,transpose
        

        if mode=='diagnostic':
            # Setup a 2D diagnostic list for the matrix and a list for the external source
            ext_source=[[] for x in range(len(self.species))]
            ret=[[[] for y in range(len(self.species))] for x in range(len(self.species))]
        elif mode in ['E','I','Sgl']:
            ret=zeros((len(self.species),len(self.species),2+3*(mode=='Sgl')))
            ext_source=zeros((len(self.species),2+3*(mode=='Sgl')))
        else:
            # Setup a matrix and vector
            ret=zeros((len(self.species),len(self.species)))
            ext_source=zeros((len(self.species)))
 
        for i in range(len(self.species)):
            ''' Walk through each row (species)'''

            for r in self.reactions:
                ''' Sort the species of each reaction into the appropriate column '''

        
                ''' Get the energy loss term '''

                Sgl=r.getS(Te,Ti,Tm,E)
                
                # TODO: what if three-particle reaction?
                bg=r.e*ne+r.p*ni # Specify density for reactions
                if mode=='R':
                    bgm=r.rate(Te,Ti,E,ne)
                    bg=r.rate(Te,Ti,E,ne)
                    ext=0
                elif mode in ['M']:#,'Sgl']:
                    bgm=(r.e*ne)*(r.p*ni)*r.rate(Te,Ti,E,ne) # Specify density for external source
                    bg=max(bg,1)*r.rate(Te,Ti,E,ne)    # Assure that auto-processes are considered
                    ext=0
                elif mode=='Sgl':
                    bgm=(r.e*ne)*(r.p*ni)*r.rate(Te,Ti,E,ne)*Sgl[:,0] # Specify density for external source
                    bg=max(bg,1)*r.rate(Te,Ti,E,ne)*Sgl[:,0]    # Assure that auto-processes are considered
                    ext=Sgl[:,1]
                elif mode=='I':
                    bgm=(r.e*ne)*(r.p*ni)*r.rate(Te,Ti,E,ne)*Sgl[-2:,0] # Specify density for external source
                    bg=max(bg,1)*r.rate(Te,Ti,E,ne)*(abs(Sgl[-2:,0])>0)    # Assure that auto-processes are considered
                    ext=0
                elif mode=='E':
                    bgm=Sgl[-2:,0] # Specify density for external source
                    bg=Sgl[-2:,0]    # Assure that auto-processes are considered
                    ext=0


                j=None # Set flag to identify external sources

                # Loop through each reaction defined in the CRM
                for rea in range(len(r.reactants)):
                    # Find the column into which the fragments goes: if background mark external
                    try:
                        j=self.slist.index(r.reactants[rea])  # Get the product species index of the correct column
                    except:
                        continue
            
                    if self.slist[i]==r.reactants[rea]:   # If the species (row index) is a reactant, r is a depletion process 
                        ''' DEPLETION '''

                        if mode=='diagnostic':
                            ''' Diagnostic matrix '''
                            ret[i][i].append('-'+str(r.r_mult[rea])+'*'+r.database+'_'+r.name+bg)   # Print the rate to the correct element

                        elif mode in ['R','M']:
                            ''' Rate coefficient matrix '''
                            ret[i,i]-=r.r_mult[rea]*bg # Calculate the Rate coefficient and store appropriately


                for frag in range(len(r.fragments)):    # Loop through the reaction fragments
                    ''' SOURCE '''

                    # Do nothing if background fragment
                 #   if r.fragments[frag] not in self.slist: 
                 #       continue
                    
                    # If fragment enters row, add in appropriate column as defined by reactant
                    if (r.fragments[frag] in self.slist) and (self.slist.index(r.fragments[frag])==i):
                            multiplier=r.f_mult[frag]**(mode not in ['Sgl','I','E'])   # Fragment multiplier
                            if j is None: # External flag triggered, store to external source
                                ''' EXTERNAL SOURCE '''

                                if mode=='diagnostic':
                                    ''' Diagnostic matrix '''
                                    ext_source[i].append('+'+str(multiplier)+'*'+r.database+'_'+r.name+bg)
                                else:
                                    try:
                                        ext_source[i]+=multiplier*bgm+ext
                                    except:
                                        ext_source[i,:]+=multiplier*bgm+ext
                                  
                            else: # No trigger of external source, store to appropriate location in matrix
                                ''' INTERNAL SOURCE '''

                                if mode=='diagnostic':
                                    ''' Diagnostic matrix '''
                                    ret[i][j].append('+'+str(multiplier)+'*'+r.database+'_'+r.name+bg)
                                else:
                                    try:
                                        ret[i,j]+=multiplier*bg+ext
                                    except:
                                        ret[i,j,:]+=multiplier*bg+ext

                                    
        return ret,ext_source



    def populate_old(self,mode,Te,ne,Ti=None,ni=None,E=0,rad=True,Sind=None,Tm=False,Ton=True,Iind=0):
        ''' Function populating a matrix according to the chosen mode 
            populate(mode,Te,ne,*keys)
            mode    -   Matrix writing mode
                            'diagnostic'    -   Creates a 2D list of reactions handles
                            'R'             -   Creates a matrix of rate coefficients (cm**3/s)
                            'M'             -   Creates a matrix of rates (s**-1)
            Te      -   Background plasma electron temperature [eV]
            ne      -   Background plasma electron density [cm**-3]

            Optional parameters
            Ti (None)   -   Background plasma ion temperature [eV]. Ti=Te assumed if None
            ni (None)   -   Background plasma ion density [cm**-3]. ni=ne assumed if None
            E (0.1)     -   Target particle energy [eV]
            rad

            Returns
            matrix,ext_source
    
            matrix      -   The matrix with the requested elements
            ext_source  -   Vector containing the external source contributions
                            to each species. Note that the external sources always are
                            returned as volumetric rates (cm**-3 s**-1), which differ from the
                            rate matrix
            
        '''
        from numpy import zeros,array,sum,transpose
        

        if mode=='diagnostic':
            # Setup a 2D diagnostic list for the matrix and a list for the external source
            ext_source=[]
            ret=[]
            for i in range(len(self.species)):
                ext_source.append([])
                ret.append([])
                for j in range(len(self.species)):
                    ret[i].append([])
        elif mode=='Sgl':
            ret=zeros((len(self.species),len(self.species),5))
            ext_source=zeros((len(self.species),5))
        elif mode in ['E','I']:
            ret=zeros((len(self.species),len(self.species),2))
            ext_source=zeros((len(self.species),2))
        else:
            # Setup a matrix and vector
            ret=zeros((len(self.species),len(self.species)))
            ext_source=zeros((len(self.species)))
 
        for i in range(len(self.species)):
            ''' Walk through each row (species)'''

            for r in self.reactions:
                ''' Sort the species of each reaction into the appropriate column '''

        
                ''' Get the energy loss term '''
                # Create array: first index Sel,SeV,Seg,Spe,SpV,Spg, second index S,ext

                Sgl=r.getS(Te,Ti,Tm,E)
                if mode in ['I','E']:
                    Sgl=Sgl[-2:]
                
                # TODO: what if three-particle reaction?
                bg=r.e*ne+r.p*ni # Specify density for reactions
                if mode!='diagnostic':
                    bgm=(r.e*ne)*(r.p*ni) # Specify density for external source
                    bg=max(bg,1)    # Assure that auto-processes are considered
                j=None # Set flag to identify external sources

                # Loop through each reaction defined in the CRM
                for rea in range(len(r.reactants)):
                    # Find the column into which the fragments goes: if background mark external
                    try:
                        j=self.slist.index(r.reactants[rea])  # Get the product species index of the correct column
                    except:
                        continue
            
                    multiplier=r.r_mult[rea] # Get the reaction multiplier for the density

                    if self.slist[i]==r.reactants[rea]:   # If the species (row index) is a reactant, r is a depletion process 
                        ''' DEPLETION '''

                        if mode=='diagnostic':
                            ''' Diagnostic matrix '''
                            ret[i][i].append('-'+str(multiplier)+'*'+r.database+'_'+r.name+bg)   # Print the rate to the correct element

                        elif mode=='R':
                            ''' Rate coefficient matrix '''
                            ret[i,i]-=multiplier*r.rate(Te,Ti,E,ne) # Calculate the Rate coefficient and store appropriately

                        elif mode=='M':
                            ''' Rate matrix '''
                            ret[i,i]-=multiplier*r.rate(Te,Ti,E,ne)*bg # Calculate the rate and store appropriately

                for frag in range(len(r.fragments)):    # Loop through the reaction fragments
                    ''' SOURCE '''
                    multiplier=r.f_mult[frag]   # Fragment multiplier

                    # Do nothing if background fragment
                    if r.fragments[frag] not in self.slist: 
                        continue
                    
                    # If fragment enters row, add in appropriate column as defined by reactant
                    elif self.slist.index(r.fragments[frag])==i:
                            if j is None: # External flag triggered, store to external source
                                ''' EXTERNAL SOURCE '''

                                if mode=='diagnostic':
                                    ''' Diagnostic matrix '''
                                    ext_source[i].append('+'+str(multiplier)+'*'+r.database+'_'+r.name+bg)

                                elif mode=='R':
                                    ''' Rate coefficient matrix '''
                                    ext_source[i]+=multiplier*r.rate(Te,Ti,E,ne)

                                elif mode=='M':
                                    ''' Rate matrix '''
                                    ext_source[i]+=multiplier*r.rate(Te,Ti,E,ne)*bgm

                                elif mode=='Sgl':
                                    ''' Energy source matrix in Greenland form '''
                                    ext_source[i,:]+=r.rate(Te,Ti,E,ne)*bgm*Sgl[:,0]+Sgl[:,1]
            
                                elif mode=='I':
                                    ''' Intensity matrix '''
                                    ext_source[i,:]+=r.rate(Te,Ti,E,ne)*bgm*Sgl[:,0]

                                elif mode=='E':
                                    ''' Intensity matrix '''
                                    ext_source[i,:]+=Sgl[:,0]


                            else: # No trigger of external source, store to appropriate location in matrix
                                ''' INTERNAL SOURCE '''

                                if mode=='diagnostic':
                                    ''' Diagnostic matrix '''
                                    ret[i][j].append('+'+str(multiplier)+'*'+r.database+'_'+r.name+bg)
                            
                                elif mode=='R':
                                    ''' Rate coefficient matrix '''
                                    ret[i,j]+=multiplier*r.rate(Te,Ti,E,ne)

                                elif mode=='M':
                                    ''' Rate matrix '''
                                    ret[i,j]+=multiplier*r.rate(Te,Ti,E,ne)*bg

                                elif mode=='Sgl': 
                                    ''' Energy source matrix in Greenland form '''
                                    ret[i,j,:]+=r.rate(Te,Ti,E,ne)*bg*Sgl[:,0]+Sgl[:,1]

                                elif mode=='I':
                                    ''' Correlation matrix '''    
                                    ret[i,j,:]+=r.rate(Te,Ti,E,ne)*bg*(abs(Sgl[:,0])>0)

                                elif mode=='E':
                                    ''' Correlation matrix '''    
                                    ret[i,j,:]+=Sgl[:,0]
                                    
        return ret,ext_source

        
    def write_matrix(self,mat,ext,char,te,ne,ti,ni,E,form='{:1.1E}'):
        ''' Writes the matrix mat to a file and outputs to stdout, if requested
            write_matrix(mat,ext,char,te,ne,ti,ni,E)
            
            mat     -   NxN matrix to be written
            ext     -   External source, N vector
            char    -   Character to be identify the matrix in the file name
            te      -   Electron temperature for header [ev]
            ne      -   Electron density for header [cm**-3]
            ti      -   Ion temperature for header [ev]
            ni      -   Ion density for header [cm**-3]

        '''
        from os import getcwd
        from datetime import datetime

        Slabels=['S_e','S_ia','S_V','S_g']

        # Create a file to write to according to char
        with open('{}/logs/{}_matrix.log'.format(self.path,char),'w') as f:  
            
            if char=='R': # Rate coefficient  matrix is being written
                f.write('''Diagnostic rate coefficient (density-independend) 
                            matrix for CRUM run in {} on {}\n
                            '''.format(getcwd(),str(datetime.now())[:-7]))
                f.write('Te={} eV, Ti={} eV, E={} eV\n'.format(te,ti,E))
            elif char=='M': # Rate matrix is being written
                f.write('''Diagnostic rate (density-dependend) matrix for CRUM 
                            run in {} on {}\n'''.format(getcwd(),
                            str(datetime.now())[:-7]))
                f.write('''Te={} eV, Ti={} eV, ne={} 1/cm**3, ni={} 
                            1/cm**3, E={} eV\n'''.format(te,ti,ne,ni,E))
            elif char=='S': # Rate matrix is being written
                f.write('''Diagnostic energy loss (density-dependend) matrix 
                            for CRUM run in {} on {}\n'''.format(getcwd(),
                            str(datetime.now())[:-7]))
                f.write('''Te={} eV, Ti={} eV, ne={} 1/cm**3, ni={} 1/cm**3, 
                            E={} eV\n'''.format(te,ti,ne,ni,E))
            
            # Create header line
            out='{}-MAT|'.format(char.upper()).rjust(10) 
            for s in self.slist:
                out+=s.rjust(10,' ')
            f.write(out+'\n'+'_'*(1+len(self.species))*10+'\n')

            # Loop through the matrix
            for l in range(len(mat)):
                if len(mat)==len(self.species):
                    # Create a tabulated file with the row species displayed
                    out=self.slist[l]+'|' 
                else:
                    # Create a tabulated file with the row species displayed
                    out=Slabels[l]+'|' 
                out=out.rjust(10,' ')
                # Add each element to the line with 10 characters reserved
                for e in mat[l,:]:
                    if e==0:
                        out+=' '*9+'_'
                    else:
                        out+=form.format(e).rjust(10,' ')
                f.write(out+'\n')
            # Write the external source to the bottom of the output
            f.write('_'*(1+len(self.species))*10+'\n')
            out='S_ext|'.rjust(10,' ')
            for s in ext:
                out+=form.format(s).rjust(10,' ')
            f.write(out)



    def DIAGNOSTIC(self):
        ''' Returns a diagnostic matrix containing lists of reactions accounted for in each element '''
        from os import getcwd
        from datetime import datetime
 
        # Get the 2D diagnostic list and external source list
        dia,ext=self.populate('diagnostic',0,'*[ne]',0,'*[ni]') 

        # Write the diagnostic matrix to the logs
        with open('{}/logs/reaction_matrix.log'.format(self.path),'w') as f:
            # Write the header line
            f.write('''Diagnostic reaction matrix for CRUM run in {} 
                        on {}\n'''.format(getcwd(),str(datetime.now())[:-7]))
            # Loop through each species
            for i in range(len(dia)):
                # Write the species header
                f.write('\n\n\n======== {} ========\n'.format(self.slist[i]))
                # Start with the sinks
                f.write('''{} DEPLETION:\n                 {}\n
                            '''.format(self.slist[i],dia[i][i]))
                # Then do the sources
                for j in range(len(dia[i])):
                    if len(dia[i][j])>0:    # Don't write anything empty
                        if i!=j:    # Don't duplicate depletion
                            f.write('''From {}:\n                 {}\n
                                        '''.format(self.slist[j],dia[i][j]))
                if len(ext[i])>0:   # Write external source last, if applicable
                    f.write('''{} EXTERNAL SOURCE:\n                 {}\n
                                        '''.format(self.slist[i],ext[i]))
        # If verbose, output the log content
        if self.verbose:
            with open('logs/reaction_matrix.log','rt') as f:
                for l in f:
                    print(l.strip())
            
            

    def R(self,Te,Ti=None,E=0.1,sparse=False,write=True):
        ''' Creates the rate coefficient matrix
            R(Te,*keys)
        
            Te              -   Background plasma electron temperature [eV]

            Optional parameters
            Ti (None)       -   Background plasma ion temperature [eV]. Ti=Te if Ti is None
            E (0.1)         -   Target particle energy [eV]
            sparse (False)  -   Switch for returning the matrix as a csc matrix
            write (True)    -   Write the rate coefficient matrix to file

            Returns
            R,ext
            R   -   Rate coefficient matrix
            ext -   External source matrix 

        '''
        from scipy.sparse import csc_matrix
        if Ti is None: Ti=Te # Check for Ti, set if necessary
        ni,ne=1,1 # Set densities to one to get rate coefficients as output
        R,ext=self.populate('R',Te,0,Ti,0,E)
    
        if write: # Write to log if requested
            self.write_matrix(R,ext,'R',Te,0,Ti,0,E)
            if self.verbose: # Print rate matrix to stdout if running verbose
                with open('logs/R_matrix.log','rt') as f:
                    for l in f:
                        print(l.strip())
        
        if sparse: R=csc_matrix(R)  # Use sparse format if requested

        return R,ext




    def M(self,Te,ne,Ti=None,ni=None,E=0.1,sparse=False,write=True):
        ''' Creates the rate  matrix
            M(Te,*keys)

            Te              -   Background plasma electron temperature [eV]

            Optional parameters
            Ti (None)       -   Background plasma ion temperature [eV]. Ti=Te if Ti is None
            E (0.1)         -   Target particle energy [eV]
            sparse (False)  -   Switch for returning the matrix as a csc matrix
            write (True)    -   Write the rate coefficient matrix to file

            Returns
            M,ext
            M   -   Rate matrix
            ext -   External source matrix 

            '''
        from scipy.sparse import csc_matrix
        if Ti is None: Ti=Te # Check for Ti, set if necessary
        if ni is None: ni=ne # Check for ni, set if necessary

        M,ext=self.populate('M',Te,ne,Ti,ni,E)

        if write:   # Write to log if requested
            self.write_matrix(M,ext,'M',Te,ne,Ti,ne,E)
            if self.verbose: # Print output if running verbose
                with open('logs/M_matrix.log','rt') as f:
                    for l in f:
                        print(l.strip())
        
        if sparse: M=csc_matrix(M) # Use sparse format if requested
    
        return M,ext
 
    def Sgl(self,Te,ne,Ti=None,ni=None,E=0.1,rad=True,Tm=False,
                        write=False,Ton=True):
        from numpy import matmul,block
        from numpy.linalg import inv
        if Ti is None: Ti=Te # Check for Ti, set if necessary
        if ni is None: ni=ne # Check for ni, set if necessary

        mat,ext=self.populate('Sgl',Te,ne,Ti,ni,E,Tm=Tm,rad=rad,Ton=Ton)    


        if write:
            title=['Sgl_el','Sgl_ia','Sgl_v','Sgl_ga','Sgl_gm']
            for i in range(5):
                self.write_matrix(mat[:,:,i],ext[:,i],title[i],
                                        Te,ne,Ti,ni,E,form='{:1.2E}')

        
        U=[ [mat[:,:,0], ext[:,0]],
            [mat[:,:,1], ext[:,1]],
            [mat[:,:,2], ext[:,2]],
            [mat[:,:,3], ext[:,3]],
            [mat[:,:,4], ext[:,4]], ]


        M,G=self.M(Te,ne,Ti,ni,E,write=write) # Get the full rate matrix

        MP=M[:self.Np,:self.Np]
        MQ=M[self.Np:,self.Np:]
        V=M[self.Np:,:self.Np]
        H=M[:self.Np,self.Np:]
        

        ret=[]
        for S in U:
            UP=S[0][:self.Np,:self.Np]
            UQ=S[0][self.Np:,self.Np:]
            UV=S[0][self.Np:,:self.Np]
            UH=S[0][:self.Np,self.Np:]

            P=UP-matmul(UH,matmul(inv(MQ),V))
            Q=UV-matmul(UQ,matmul(inv(MQ),V))
            
            ret.append([block( [ [P],[Q]]),S[1]])

            
        return ret

    def S(self,Te,ne,Ti=None,ni=None,E=0.1,rad=True,
                                Tm=False,write=False,Ton=True):
        if Ti is None: Ti=Te # Check for Ti, set if necessary
        if ni is None: ni=ne # Check for ni, set if necessary

        mat,ext=self.populate('Sgl',Te,ne,Ti,ni,E,Tm=Tm,rad=rad,Ton=Ton)    


        if write:
            title=['Sgl_el','Sgl_ia','Sgl_v','Sgl_ga','Sgl_gm']
            for i in range(5):
                self.write_matrix(mat[:,:,i],ext[:,i],title[i],
                                    Te,ne,Ti,ni,E,form='{:1.2E}')

        return  [   [mat[:,:,0], ext[:,0]],
                    [mat[:,:,1], ext[:,1]],
                    [mat[:,:,2], ext[:,2]],
                    [mat[:,:,3], ext[:,3]],
                    [mat[:,:,4], ext[:,4]], ]


    def I(self,Te,ne,Ti=None,ni=None,E=0.1,
                                rad=True,Tm=False,write=False,Ton=True):
        ''' Creates a radiaiton intensity matrix '''
        if Ti is None: Ti=Te # Check for Ti, set if necessary
        if ni is None: ni=ne # Check for ni, set if necessary

        mat,ext=self.populate('I',Te,ne,Ti,ni,E,Tm=Tm,rad=rad,Ton=Ton)    
        if write:
            title=['Ia','Im']
            for i in range(5):
                self.write_matrix(mat[:,:,i],ext[:,i],title[i],
                                            Te,ne,Ti,ni,E,form='{:1.2E}')


        return  [   [mat[:,:,0], ext[:,0]],
                    [mat[:,:,1], ext[:,1]] ]


    def E(self,Te,ne,Ti=None,ni=None,E=0.1,rad=True,
                                Tm=False,write=False,Ton=True):
        ''' Creates a radiaiton intensity matrix '''
        if Ti is None: Ti=Te # Check for Ti, set if necessary
        if ni is None: ni=ne # Check for ni, set if necessary

        mat,ext=self.populate('E',Te,ne,Ti,ni,E,Tm=Tm,rad=rad,Ton=Ton)    


        if write:
            title=['Ea','Em']
            for i in range(5):
                self.write_matrix(mat[:,:,i],ext[:,i],
                                title[i],Te,ne,Ti,ni,E,form='{:1.2E}')


        return  [   [mat[:,:,0], ext[:,0]],
                    [mat[:,:,1], ext[:,1]] ]



    def full_nt(self,Te,ne,t,gl,Ti=None,ni=None,E=0.1,n=None,Sext=True):
        ''' Solves the full NxN problem 
            full_nt(Te,ne,t,*keys)
            
            Te      -   Background plasma electron temperature [eV]
            ne      -   Background plasma electron density [cm**-3]

            Optional parameters
            Ti (None)   -   Background plasma ion temperature [eV]. Ti=Te assumed if None
            ni (None)   -   Background plasma ion density [cm**-3]. ni=ne assumed if None
            E (0.1)     -   Target particle energy [eV]
            Sext (True) -   Include external source (from background plasma reactions into CRM species)
            n (None)    -   Initial distribution of particeles, taken as n0 specified in input if None

            Returns
            ivp_solve bunch object containing the time-dependent ODE system solution
        '''
        from scipy.integrate import solve_ivp 

        # Use input n0 as default unless explict n0 requested
        if n is None: n=self.n0() 

        if gl is True:
            # Set up Greenland model 
            mat,ext,nP0p=self.gl_crm(*self.M(Te,ne,Ti,ni,E,write=False),Sext,n)
            n=nP0p
        else:
            mat,ext=self.M(Te,ne,Ti,ni,E,write=False) # Get the full rate matrix

        ext=(Sext is True)*ext  # Set source strength
        # Solve and return
        return solve_ivp(lambda x,y: self.ddt(x,y,mat,ext),
                                            (0,t),n,method='LSODA')





        
    def dEdt(self,t,Te,ne,Ti=None,ni=None,E=0.1,Tm=False,rad=True,
                Sext=True,write=False,gl=True,n=None,Qres=True,Ton=True):
        from numpy import block,zeros,matmul,reshape,sum
        from  numpy.linalg import inv
        from scipy.integrate import solve_ivp


        # TODO: Combine with full_nt  
        N=len(self.species)
        Np=self.Np
        Nq=N-Np  


        

        M,G=self.M(Te,ne,Ti,ni,E,write=write) # Get the full rate matrix

        ret=[]
        if gl is True:
            if n is None:
                n0=self.n0()
            else:
                n0=n
            # Create block matrix from M
            MP=M[:self.Np,:self.Np]
            MQ=M[self.Np:,self.Np:]
            V=M[self.Np:,:self.Np]
            H=M[:self.Np,self.Np:]

            Meff,GPp,nP0p=self.gl_crm(M,G,Sext,n0)
        
            U=self.Sgl(Te,ne,Ti,ni,E,rad,Tm,write=write,Ton=Ton)


            if Qres is True:

                mat=block(  [   [zeros((N,5*N)),     U[0][0]],
                                [zeros((N,5*N)),     U[1][0]],
                                [zeros((N,5*N)),     U[2][0]],
                                [zeros((N,5*N)),     U[3][0]],
                                [zeros((N,5*N)),     U[4][0]],
                                [zeros((Np,5*N)),     Meff]   ])
                if n is None:
                    n=zeros((len(mat),))
                    n[-Np:]=nP0p
                ext=block([U[0][1], U[1][1], U[2][1], U[3][1], U[4][1], GPp])

            else:

                mat=block(  [   [zeros((5,)),      sum(U[0][0],axis=0)],
                                [zeros((5,)),      sum(U[1][0],axis=0)],
                                [zeros((5,)),      sum(U[2][0],axis=0)],
                                [zeros((5,)),      sum(U[3][0],axis=0)],
                                [zeros((5,)),      sum(U[4][0],axis=0)],
                                [zeros((Np,5)),   Meff]   ])

                if n is None:
                    n=zeros((len(mat),))
                    n[-Np:]=nP0p
                ext=block([ sum(U[0][1],axis=0), sum(U[1][1],axis=0), 
                            sum(U[2][1],axis=0), sum(U[3][1],axis=0), 
                            sum(U[4][1],axis=0), GPp                    ])
                

     

        else:
             
            if Qres is True:
                U=self.S(Te,ne,Ti,ni,E,rad,Tm,write=True,Ton=Ton)
    

                mat=block(  [   [zeros((N,5*N)),        U[0][0]],
                                [zeros((N,5*N)),        U[1][0]],
                                [zeros((N,5*N)),        U[2][0]],
                                [zeros((N,5*N)),        U[3][0]],
                                [zeros((N,5*N)),        U[4][0]],
                                [zeros((N,5*N)),        M]   ])
                if n is None:
                    n=zeros((len(mat),))
                    n[-N:]=self.n0()
                ext=block([U[0][1], U[1][1], U[2][1], U[3][1], U[4][1], G])

            else:
                U=self.S(Te,ne,Ti,ni,E,rad,Tm,write=True,Ton=Ton)

                mat=block(  [   [zeros((5,)),      sum(U[0][0],axis=0)],
                                [zeros((5,)),      sum(U[1][0],axis=0)],
                                [zeros((5,)),      sum(U[2][0],axis=0)],
                                [zeros((5,)),      sum(U[3][0],axis=0)],
                                [zeros((5,)),      sum(U[4][0],axis=0)],
                                [zeros((N,5)),    M]   ])

                if n is None:
                    n=block([   zeros((5,)), self.n0()    ] )

                n=n.reshape((N+5,))

                ext=block([ sum(U[0][1],axis=0), sum(U[1][1],axis=0),
                            sum(U[2][1],axis=0), sum(U[3][1],axis=0), 
                            sum(U[4][1],axis=0), G                      ])

        return solve_ivp(lambda x,y: self.ddt(x,y,mat,ext),(0,t),
                                            n,'LSODA',dense_output=True)

        

    #def gl_E(self,mat,ext,Sext=True,n=None,matrices=False):


    def gl_crm(self,mat,ext,Sext=True,n=None,matrices=False):
        ''' Returns the P-space matrices according to Greenland 2001
            gl_crm(ne,*keys)

            Optional parameters
            Sext (True) -   Include external source (from background plasma reactions into CRM species)
            n (None)    -   Initial distribution of particeles, taken as n0 specified in input if None
            matrices (False)    -   Switch determining whether to return the CRM or the matrices

            Returns (matrices=False)
            Meff,GPp,nP0p
            Meff        -   Effective rate matrix, Np x Np
            GPp         -   The modified external source, accounting for 
                            reactions turning Q-species into P-species
            nP0p        -   The modified initial density for the P-space, 
                            accounting for Q->P space reactions

            Returns (matrices=True)
            M           -   The full rate matrix
            T           -   The normalized, right eigenvector matrix
            D           -   The diagonalized eigenvalue matrix

        '''
        from numpy import matmul,diag,real
        from numpy.linalg import inv,eig

        # Use input n0 as default unless explict n0 requested
        if n is None: n=self.n0() 
        

        # Create block matric from M
        MP=mat[:self.Np,:self.Np]
        MQ=mat[self.Np:,self.Np:]
        V=mat[self.Np:,:self.Np]
        H=mat[:self.Np,self.Np:]

        # Calculate Meff
        Meff=(MP-matmul(matmul(H,inv(MQ)),V))

        # Diagonalize M
        eigs,T=eig(mat)

        if matrices is True:
            return mat,T,real(eigs)
        # Order the eigenvalues and vectors in increasing magnitude
        # TODO should this be done for all the eigenvectors or separately   
        # between the Q and P space??
        # --> Does not seem to affect the solution
        eigind=abs(eigs[self.Np:]).argsort()[::1] 
        e1=eigs[self.Np:]
        T1=T[:,self.Np:]
        e1=e1[eigind]
        T1=T1[:,eigind]

        #eigs=eigs[eigind]
        #T=T[:,eigind]
        eigs[self.Np:]=e1
        T[:,self.Np:]=T1
        # Only use real parts TODO is this OK??
        eigs=real(eigs)
        D=diag(eigs)

        # Create block matrix
        TQ=T[self.Np:,self.Np:]
        TP=T[:self.Np,:self.Np]
        delta=T[self.Np:,:self.Np]
        Delta=T[:self.Np,self.Np:]
        

        # Calculate P-space CRM
        GPp=real((Sext is True)*ext[:self.Np]-matmul(matmul(Delta,inv(TQ)),ext[self.Np:]))
        nP0p=real(n[:self.Np]-matmul(matmul(Delta,inv(TQ)),n[self.Np:]))
        
        if matrices is False:
            return Meff,GPp,nP0p

  
    def ddt(self,t,n,mat,ext):
        ''' Returns the time-derivative of the density for the CRM i
            dndt(t,n,mat,ext)
        
            t   -   Time of evaluation [s]
            n   -   Vector of initial density distribution [cm**-3]
            mat -   Rate matrix
            ext -   External source vector
        '''
        from numpy import matmul

        return matmul(mat,n)+ext

    def evaluate_CRM(self,Te,ne,Ti=None,ni=None,E=0.1,printout=True):
        ''' Evaluates the current CRM '''
        from numpy import matmul
        from numpy.linalg import inv
        # Get matrices and eigenvalues/-vectors
        M,T,D=self.gl_crm(*self.M(Te,ne,Ti,ni,E,write=False),matrices=True) 
        TQ,delta=T[self.Np:,self.Np:],T[self.Np:,:self.Np]
        tau=1/abs(D)
        tauPmin=min(tau[:self.Np])
        tauQmax=max(tau[self.Np:])
        norm=matmul(inv(TQ),delta)
        maxnorm=max([sum(abs(norm[:,i])) for i in range(norm.shape[1])])
        maxdelta=max([sum(abs(delta[:,i])) for i in range(delta.shape[1])])
        if printout is True:
            print("Validity of CRM with Np={}".format(self.slist[:self.Np]))
            print("    max(|delta|)={:.2E}".format(maxdelta))
            print("    max(|inv(T_Q)*delta|)={:.2E}".format(maxnorm))
            print("    min(tau_P)={:.2E}".format(tauPmin))
            print("    max(tau_Q)={:.2E}".format(tauQmax))
        else:
            return [maxdelta,maxnorm,tauPmin,tauQmax]

    
    def generate_CRM(self,Te,ne,kappa,Ti=None,ni=None,E=0.1,n=None,
                                                    Sext=True,epsilon=1):
        ''' Generates the optimal CRMs per Greenland 2001 '''
        from numpy import array,where,matmul
        from numpy.linalg import inv
        # Get matrices and eigenvalues/-vectors
        M,T,D=self.gl_crm(*self.M(Te,ne,Ti,ni,E,write=False),
                                        Sext=Sext,n=n,matrices=True) 
        # Construct indicator matrix
        I=abs(T)
        I[I<=kappa]=0
        I[I>kappa]=1

        # Create Z-array
        Z=[]
        for i in I: # Loop over species
            N=0
            while i[N]==0: # Count leading zeros
                N+=1
                # Break if all species are 0
                if N==len(i):
                    break
            Z.append(N)

        # Get increasing order of leading zeros
        Z=array(Z)
        Zind=Z.argsort()[::1] 
        # Check whether indicator matrix is singular
        if max(Z)==M.shape[1]:
            print('WARNING! Matrix I is singular')

        # Identify a boolean array for valid Np
        Np=[Z[Zind][i]>=i for i in range(len(Z))]


        di=range(len(Z))-Zind # Track change in index space
        sorted_D=D[Zind]
        sorted_species=[self.slist[i] for i in Zind]
        rej=[]
        for i in range(len(Z)):
            if Np[i]==True:
                if i+1==len(Z): continue # Full model, uninteresting
                TQ,delta=T[i+1:,i+1:],T[i+1:,:i+1]
                norm=matmul(inv(TQ),delta)
                maxnorm=max([sum(abs(norm[:,i])) for i in range(norm.shape[1])])
                if maxnorm<epsilon:
                    print('''Possible CRM NP={} model with TQinv*delta={:.2E}, 
                            tauP={:.2E} s, tauQ={:.2E} s, and  P-space {} 
                            found.'''.format(i+1,maxnorm,
                            1/abs(min(sorted_D[:i+1])),
                            1/abs(max(sorted_D[i+1:])),sorted_species[:i+1])) 
                else:
                    rej.append(maxnorm)
        if len(rej)>0:
            print('Minimum rejected maxnorm={:.2E}'.format(min(rej)))



        
    def intensity(self,Te,ne,Ti=None,ni=None,E=0.1,units='v',norm=True,
                            write=False,Sext=True,n=None,reassociation=1e3):
        ''' Returns a list of 2x(NxN) vector with spectroscopic units/intensity of [Ia, Im]'''
        from numpy import reshape,zeros,multiply,where
        from numpy import transpose,matmul,array,mean,all,diag
        from scipy.optimize import minimize
        from scipy.integrate import solve_ivp
        from matplotlib.pyplot import figure


        '''
        def dt(t,n,mat,ext):
            m=matmul(mat,n)
            ext[0]=-m[0]
            return matmul(mat,n)+ext



        if n is None: n=self.n0() # Use input n0 as default unless explict n0 requested
        NN=len(self.species)*len(self.species)

        mat,ext=self.M(Te,ne,Ti,ni,E,write=False) # Get the full rate matrix
       
        # Supply molecules to achieve SS
        ext[1]=sum(n)/2
        # Simulate to SS
        ss_sol=solve_ivp(lambda x,y: dt(x,y,mat,ext),(0,10),n,'LSODA',dense_output=True)

        ss_sol=diag(ss_sol.y[:,-1])
        '''

        def dt(t,n,mat,ext,n0):
            '''
            ntot=self.totparticles(n)
            dn=self.totparticles(n)-self.totparticles(n0)
            ext[1]=-dn*fac
            '''
            #fac=0*1e3#-1e-2
            '''ext[0]=-(matmul(mat,n))[0]'''
            
            #ext[1]=-0.5*ext[0]
            #ext[1]=-ext[0]*0.1
            #ext[0]=-fac*n[0]
            #ext[1]=-0.99*self.totpart(ext)
            #ext[1]=-0.5*ext[0]
            return matmul(mat,n)+ext

        def ss(t,Y):
            global ss_n 
            if norm(Y-ss_n)<7e5:
                return 0
            else:
                return 1
        ss.terminal=True

        def findfac(fac,n,mat,ext,n0):
            nend_sum=self.totparticles(solve_ivp(lambda x,y: dt(x,y,mat,ext,n,
                            fac),(0,1),n,'LSODA',dense_output=True).y[:,-1])
            
            return abs(self.totparticles(n0)-nend_sum)
            
            


        if n is None: n=self.n0() # Use input n0 as default unless explict n0 requested
        NN=len(self.species)*len(self.species)

        mat,ext=self.M(Te,ne,Ti,ni,E,write=False) # Get the full rate matrix
       
        ''' Assume re-association to achieve SS 
        mat[0,0]=-reassociation
        mat[1,0]=-mat[0,0]*0.5 # Conserve nucleii
    
        fac=minimize(findfac,5e4,(n,mat,ext,n))
        '''

        '''
        psor=-4.119965e20 
        psorgc=3.26090e20
        n[0]=4.846e13
        n[1]=2.925e13
        ext[0]=psorgc*1e-6-8.552e14
        ext[1]=6.33578e14   
        mat[0,0]=(psor*1e-6)/ne
        '''
        # Supply molecules to achieve SS
        '''ext[1]=sum(n)/0.2'''
        # Simulate to SS
        ss_sol=diag(solve_ivp(lambda x,y: dt(x,y,mat,ext,n),(0,1),
                                    n,'LSODA',dense_output=True).y[:,-1])
        #return







        #def ss(n,mat):
        #    return matmul(mat,n)

        #if n is None: n=self.n0() # Use input n0 as default unless explict n0 requested
        #NN=len(self.species)*len(self.species)

        #mat,ext=self.M(Te,ne,Ti,ni,E,write=write) # Get the full rate matrix
        # Take 2 'time-steps' for initial guess
        #for i in range(2):
        #    n=n+(matmul(mat,n)+ext*(Sext is True))*1e-9


        x=self.E(Te,ne,Ti,Te,E,write=False)
        y=self.I(Te,ne,Ti,Te,E,write=False)
   

        if write:
            self.write_matrix(y[0][0],n,'Ia',Te,ne,Ti,ni,E,form='{:1.2E}')
            self.write_matrix(y[1][0],n,'Im',Te,ne,Ti,ni,E,form='{:1.2E}')

        
        ret=[]
        #ss_sol=diag((fsolve(ss,n,(mat))))
        for i in range(2):
            arr=zeros((NN,2))

            # Use the steady-state solution to find transitions 
            Y=matmul(y[i][0],ss_sol)
            arr[:,0]=reshape(x[i][0],(1,NN))
            arr[:,1]=reshape(Y,(1,NN))

            # Drop non-radiative entries
            arr=transpose(arr[~all(abs(arr)< 1e-8,axis = 1)])
        
            if units=='ev':     arr[0,:]=1*arr[0,:]
            elif units=='v':    arr[0,:]=1e7/(1239.84193/arr[0,:])
            elif units=='f':    arr[0,:]=241.798*arr[0,:]
            elif units=='l':    arr[0,:]=1239.84193/arr[0,:]
            elif units=='Å':    arr[0,:]=12398.4193/arr[0,:]


            if norm is True:
                nrm=sum(arr[1,:])

            else:
                nrm=1
            arr[1,:]=arr[1,:]/nrm

            ret.append(arr)

        return ret
            

    def species_rate(self,key,Te,ne,Ti=None,ni=None,E=0.1):
        if Ti is None: Ti=Te
        if ni is None: ni=ne

        ret=0
        for r in self.reactions:
            mult=1
            if key in r.reactants:
                if r.type=='COEFFICIENT':
                    mult=1/ne
                ret+=r.rate(Te,Ti,E,ne)*mult

        return ret








    def n0(self):
        from numpy import zeros
        ret=zeros((len(self.species),))
        for s in range(len(self.species)):
            try:
                ret[s]=self.species[self.slist[s]]['n'] 
            except:
                pass
        return ret


    def get_reaction(self,database,name):
        for r in self.reactions:
            if name==r.name:
                if r.database==database:
                    return r        
        return None

    def get_rate(self,database,name,T,n,E=0.1):
        return self.get_reaction(database,name).rate(T,T,E,n)




    def totparticles(self,arr,V=1):
        ''' Returns the total number of nuclei for each time-step or of the full vector
            totparticles(arr,*keys)
            
            arr     -   1D or 2D array contining the densities of CRM species
                        The first dimension are the densities defined according
                        to the species array [cm**-3]
                        The second dimension should be the time-steps

            Optional parameters
            V (1)   -   Volume of the test box [cm**3]

            Function assumes any molecular species to contain the string 'H2' in 
            their identifier!
        '''
        from numpy import zeros

        
        try:    # If matrix 
            ret=zeros((arr.shape[1],))
            for j in range(arr.shape[1]):   # Loop through eac time-step
                for i in range(len(arr)):   # Loop through each species in arr
                    ret[j]+=(1+('H2' in self.slist[i]))*arr[i,j]
        except: # If vector
            ret=0
            for i in range(len(arr)): # Loop through each species in arr
                ret+=(1+('H2' in self.slist[i]))*arr[i]

        return ret

    def totmol(self,arr):
        from numpy import zeros

        
        try:    # If matrix 
            ret=zeros((arr.shape[1],))
            for j in range(arr.shape[1]):   # Loop through eac time-step
                for i in range(len(arr)):   # Loop through each species in arr
                    ret[j]+=(('H2' in self.slist[i]))*arr[i,j]
        except: # If vector
            ret=0
            for i in range(len(arr)): # Loop through each species in arr
                ret+=(('H2' in self.slist[i]))*arr[i]
        return ret



















