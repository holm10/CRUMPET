'''=====================================================================
========================================================================
====================COLLISIONAL-RADIATIVE MODEL=========================
========================================================================
====================================================================='''





'''=====================================================================
========================================================================
====================COLLISIONAL-RADIATIVE MODEL=========================
========================================================================
====================================================================='''

class CRM:
    # Collisional-radiative model
    # Contains:
    #   PARAMETERS
    #   FUNCTIONS:
    #   

    def __init__(self, path='.',fname='input.dat'):
#,amjuel='amjuel.tex',hydhel='hydhel.tex',h2vibr='h2vibr.tex',ADAS='ic#h0-1.dat',UE='ehr1.dat'):
        import read_rate_data as rrd
        from importlib import reload
        from numpy import zeros
        from os import mkdir,getcwd
        from datetime import datetime
        reload(rrd)
        self.verbose=False

        def file2list(path,fname):
            ''' Returns a list of lines based on fname with comment and empty lines removed as well as the card and subcard locations in the list'''
            ret=[]
            # Parse out comments and empty lines and store data to memory
            with open('{}/{}'.format(path,fname)) as f:  # Open the file
                for l in f:
                    if len(l.split())==0:
                        continue
                    elif l.strip()[0]=='#':
                        continue
                    else:
                        ret.append(l)
        
            # Locate the cards
            cards=[i for i,x in enumerate(ret) if x[0:2]=='**']
            cards.append(len(ret))
            # Locate the subcards
            subcards=[i for i,x in enumerate(ret) if x[0]=='*']
            subcards.append(len(ret))

            return ret,cards,subcards

        def XY2num(string,X=None,Y=None):
            ''' Replaces 'X' or 'Y' with X and Y, respectively - or both, depending on which is defined '''
            X=str(X)
            Y=str(Y)
            if Y is None:
                return string.replace('X',X)
            elif X is None:
                return string.replace('Y',Y)
            else:
                return string.replace('X',X).replace('Y',Y)

        def rf(string):
            ''' Split reaction into reactants and fragments '''
            return string.split(' > ')

        # Read the input file into a list
        lines,cards,subcards=file2list(path,fname)


        ''' LOOP THROUGH THE DEFINED CARDS '''
        for i in range(len(cards)-1):
            if lines[cards[i]].split()[1].upper()=='SPECIES':
                ''' Store species from species card to object '''
                self.species=[]
                for j in range(cards[i]+1,cards[i+1]):
                    self.species.append(lines[j].split()[0].strip())
                # Set up array for initial densities
                self.n0=zeros((len(self.species),))

            elif lines[cards[i]].split()[1].upper()=='REACTIONS':
                '''  Store the reactions requested to temporary reaction list '''
                reactions=[]
                for j in range(cards[i]+1,cards[i+1],2):
                    if lines[j].split()[1].upper()=='CUSTOM':
                        reactions.append( [lines[j].split()[1], lines[j+1]])
                    else:
                        reactions.append([  lines[j].split()[1],
                                            *[x.split(' + ') for x in rf(lines[j+1])],
                                         ])
            elif lines[cards[i]].split()[1].upper()=='RATES':
                ''' Store the locations of the standard rate data files '''
                for j in range(cards[i]+1,cards[i+1]):
                    [db,nm]=lines[j].split()
                    db=db.upper()
                    if db=='ADAS':
                        ADAS=nm
                    elif db=='UE':
                        UE=nm
                    elif db=='AMJUEL':
                        amjuel=nm
                    elif db=='HYDHEL':
                        hydhel=nm
                    elif db=='H2VIBR':
                        h2vibr=nm
                    else:
                        print('Unrecognized database type {}. Ignoring.'.format(db))
                        continue

            elif lines[cards[i]].split()[1].upper()=='SETTINGS':
                ''' Setup the CRM '''
                for j in range(cards[i]+1,cards[i+1]):
                    try:
                        [c,setting,value]=lines[j].split()
                    except:
                        [c,setting]=lines[j].split()
                        
                    setting=setting.upper()
                    if c!='*':
                        continue
                    elif setting=='VMAX':
                        vmax=int(value)
                    elif setting=='NMAX':
                        nmax=int(value)
                    elif setting=='VERBOSE':
                        self.verbose=bool(int(value))
                    elif setting=='NP':
                        self.Np=int(value)
                    elif setting=='N0':
                        for k in range(j+1,subcards[subcards.index(j)+1]):
                            self.n0[self.species.index(lines[k].split()[0])]=float(lines[k].split()[1])
                                                

                    else:
                        print('Unrecognized setting {}. Ignoring.'.format(setting))
                        continue

            else:
                ''' Unknown card, abort '''
                print('Unknown card "{}": aborting!'.format(lines[cards[i]].split()[1]))
                return

        ''' READ EIRENE DATA '''
        ratedata=rrd.RATE_DATA(amjuel,hydhel,h2vibr,ADAS,UE,path)


        # Create list of reaction objects
        self.reactions=[]
        for r in reactions:
            # Split the definitions into names and databases
            database=r[0].split('_')[0]
            try:
                ID=r[0].split('_')[1]
            except:
                pass

            # Check whether we have vibrationally resolved data for molecules
            if 'XvY' in ID:
                origID=ID
                for x in range(vmax+1):
                    for y in [-1,1]:
                        if x+y in range(vmax+1):
                            ID=XY2num(origID,x,x+y)
                            self.reactions.append(REACTION(    ID,
                                                                    database,
                                                                    [XY2num(i,x,x+y) for i in r[1]],
                                                                    [XY2num(i,x,x+y) for i in r[2]],
                                                                    ratedata.get_coeff(database,ID),
                                                                    'RATE'
                                                                    
                                                                ))
            elif 'Xl' in ID:
                origID=ID
                for x in range(vmax+1):
                    ID=XY2num(origID,x)
                    self.reactions.append(REACTION(    ID,
                                                            database,
                                                            [XY2num(i,x) for i in r[1]],
                                                            [XY2num(i,x) for i in r[2]],
                                                            ratedata.get_coeff(database,ID),
                                                            'RATE'
                                                        ))


            elif database.upper()=='CUSTOM':
                # Store data in data
                data,_,subcards=file2list(path,r[1].strip())
                database_use=data[0].strip()    # Use the first line as database

                # Loop through the cards
                for i in range(len(subcards)-1):
                    subc=subcards[i] # Helper
                    fit=data[subc].split()[1].upper()
                    name=data[subc].split()[2]
                    reactants,fragments=rf(data[subc+1])
                    # Identify the type of reaction we are working with
                    if fit=='RATE':
                        # Vibrationally dependent rate
                        if 'X' in name:
                            for j in range(vmax+1):
                                self.reactions.append(REACTION(    XY2num(name,j),
                                                                        database_use,
                                                                        XY2num(reactants,j).split(' + '),
                                                                        XY2num(fragments,j).split(' + '),
                                                                        [float(x) for x in data[subc+3+j*2].split()],
                                                                        fit
                                                                    )
                                                     )


                    elif fit=='COEFFICIENT':
                        if 'X' in name:
                            for j in range(vmax+1):
                                self.reactions.append(REACTION(XY2num(name,j),
                                                                    database_use,
                                                                    XY2num(reactants,j).split(' + '),
                                                                    XY2num(fragments,j).split(' + '),
                                                                    float(data[subc+2+j].split()[1]),
                                                                    fit
                                                                   )
                                                     )
                        else:
                            self.reactions.append(REACTION(    name,
                                                                    database_use,
                                                                    reactants.split(' + '),
                                                                    fragments.split(' + '),
                                                                    float(data[subc+2]),
                                                                    fit
                                                                )
                                                 )
                    elif fit=='SIGMA':
                        self.reactions.append(REACTION(    name,
                                                                database_use,
                                                                reactants.split(' + '),
                                                                fragments.split(' + '),
                                                                [float(x) for x in data[subc+2].split()],
                                                                fit
                                                            )
                                             )
                    else:
                        print('Reaction type "{}" not recognized! Aborting.'.format(data[subcards[i]].split()[1]))
                        return
                    

            elif database.upper()=='ADAS':
                for x in range(1,nmax+1):
                    if ID.upper()=='EXCITATION':
                        rn=range(x+1,nmax+1)
                        fit='ADAS'
                        Tarr=ratedata.get_coeff(database,'T')
                    elif ID.upper()=='RELAXATION':
                        rn=range(1,x)
                        fit='COEFFICIENT'
                        Tarr=0
                    for y in rn:
                            self.reactions.append(REACTION(    '{}_{}-{}'.format(ID,x,y),
                                                                    database,
                                                                    [XY2num(i,x,y) for i in r[1]],
                                                                    [XY2num(i,x,y) for i in r[2]],
                                                                    ratedata.get_coeff(database,'{}-{}'.format(x,y)),
                                                                    fit,
                                                                    Tarr
                                                                )
                                                 )

            elif database.upper()=='UE':
                self.reactions.append(REACTION(    ID,
                                                        database,
                                                        r[1],
                                                        r[2],
                                                        ratedata.get_coeff(database,ID),
                                                        'UE'
                                                   )
                                     )
                            
                    

            elif database.upper() in ['HYDHEL','AMJUEL','H2VIBR']:
                self.reactions.append(REACTION(    ID,
                                                        database,
                                                        r[1],
                                                        r[2],
                                                        ratedata.get_coeff(database,ID),
                                                        'RATE'
                                                    )
                                     )
            
            else:
                print('Database "{}" not recognized! Aborting.'.format(database))
                return

        try:
            mkdir('logs')
        except:
            pass


        # Write logs and output if requested
        with open('logs/setup.log','w') as f:
            f.write('CRUM run in {} on {}\n'.format(getcwd(),str(datetime.now())[:-7]))
            f.write('Defined species:\n')
            for i in self.species:
                f.write('    {}\n'.format(i))

            f.write('Defined reactions:\n')
            for r in self.reactions:
                f.write('{}\n'.format(r.print_reaction()))

        if self.verbose:
            with open('logs/setup.log','rt') as f:
                for l in f:
                    print(l.strip())

        self.DIAGNOSTIC()

        

    def populate(self,mode,Te,ne,Ti=None,ni=None,E=0):
        from numpy import zeros
        ''' Function populating a matrix according to the chosen mode '''
        if mode=='diagnostic':
            ext_source=[]
            ret=[]
            for i in range(len(self.species)):
                ext_source.append([])
                ret.append([])
                for j in range(len(self.species)):
                    ret[i].append([])
        else:
            ext_source=zeros((len(self.species),))
            ret=zeros((len(self.species),len(self.species)))
 
        for i in range(len(self.species)):
            ''' Walk through each row '''

            for r in self.reactions:
                ''' Sort the species of each reaction into the appropriate column '''
                
                # TODO: what if three-particle reaction?
                bg=('e' in r.reactants)*ne+('p' in r.reactants)*ni # Specify background
                if mode!='diagnostic':
                    bgm=(('e' in r.reactants)*ne)*(('p' in r.reactants)*ni) # Specify background
                    bg=max(bg,1)    # Assure that auto-processes are considered
                j=None

                # Reactants for depletion/external
                for rea in range(len(r.reactants)):
                    # Find the column into which the fragments goes: if background mark external
                    try:
                        j=self.species.index(r.reactants[rea])
                    except:
                        continue
                    multiplier=r.r_mult[rea]

                    if self.species[i]==r.reactants[rea]:
                        ''' DEPLETION '''

                        if mode=='diagnostic':
                            ''' Diagnostic matrix '''
                            ret[i][i].append('-'+str(multiplier)+'*'+r.database+'_'+r.name+bg)   # Print the rate to the correct element

                        elif mode=='R':
                            ''' Rate coefficient matrix '''
                            ret[i,i]-=multiplier*r.rate(Te,Ti,E,ne)

                        elif mode=='M':
                            ''' Rate matrix '''
                            ret[i,i]-=multiplier*r.rate(Te,Ti,E,ne)*bg

                for frag in range(len(r.fragments)):
                    ''' SOURCE '''
                    multiplier=r.f_mult[frag]   # Fragment multiplier

                    # Do nothing if background fragment
                    if r.fragments[frag] not in self.species: 
                        continue
                    
                    # If fragment enters row, add in appropriate column as defined by reactant
                    elif self.species.index(r.fragments[frag])==i:
                            if j is None:
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

                            else:
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

        return ret,ext_source

        
    def write_matrix(self,mat,ext,char,te,ne,ti,ni,E):
        ''' Outputs the numerical matrix mat to stdout '''
        from os import getcwd
        from datetime import datetime

        with open('logs/{}_matrix.log'.format(char),'w') as f:
            if char=='R':
                f.write('Diagnostic rate coefficient (density-independend) matrix for CRUM run in {} on {}\n'.format(getcwd(),str(datetime.now())[:-7]))
                f.write('Te={} eV, Ti={} eV, E={} eV\n'.format(te,ti,E))
            elif char=='M':
                f.write('Diagnostic rate (density-dependend) matrix for CRUM run in {} on {}\n'.format(getcwd(),str(datetime.now())[:-7]))
                f.write('Te={} eV, Ti={} eV, ne={} 1/cm**3, ni={} 1/cm**3, E={} eV\n'.format(te,ti,ne,ni,E))
            

            out='{}-MATRIX |'.format(char.upper())
            for s in self.species:
                out+=s.rjust(10,' ')
            f.write(out+'\n'+'_'*(1+len(self.species))*10+'\n')
            for l in range(len(mat)):
                out=self.species[l]+'|'
                out=out.rjust(10,' ')
                for e in mat[l,:]:
                    if e==0:
                        out+=' '*9+'_'
                    else:
                        out+='{:1.1E}'.format(e).rjust(10,' ')
                f.write(out+'\n')
            f.write('_'*(1+len(self.species))*10+'\n')
            out='S_ext|'.rjust(10,' ')
            for s in ext:
                out+='{:1.1E}'.format(s).rjust(10,' ')
            f.write(out)



    def DIAGNOSTIC(self):
        ''' Returns a diagnostic matrix containing lists of reactions accounted for in each element '''
        from os import getcwd
        from datetime import datetime
 
        dia,ext=self.populate('diagnostic',0,'*[ne]',0,'*[ni]')

        with open('logs/reaction_matrix.log','w') as f:
            f.write('Diagnostic reaction matrix for CRUM run in {} on {}\n'.format(getcwd(),str(datetime.now())[:-7]))
            for i in range(len(dia)):
                f.write('\n\n\n======== {} ========\n'.format(self.species[i]))
                f.write('{} DEPLETION:\n                 {}\n'.format(self.species[i],dia[i][i]))
                for j in range(len(dia[i])):
                    if len(dia[i][j])>0:
                        if i!=j:
                            f.write('From {}:\n                 {}\n'.format(self.species[j],dia[i][j]))
                if len(ext[i])>0:
                    f.write('{} EXTERNAL SOURCE:\n                 {}\n'.format(self.species[i],ext[i]))

        if self.verbose:
            with open('logs/reaction_matrix.log','rt') as f:
                for l in f:
                    print(l.strip())
            
            

    def R(self,Te,Ti=None,E=0.1,sparse=False,write=True):
        from scipy.sparse import csc_matrix
        # Returns the rate coefficient matrix at temperature Te and density ne
        # Allow for different ion and electron parameters: if not defined,set to equal
        if Ti is None:
            Ti=Te
        ni,ne=1,1 
        R,ext=self.populate('R',Te,0,Ti,0,E)
    
        if write:
            self.write_matrix(R,ext,'R',Te,0,Ti,0,E)
            if self.verbose:
                with open('logs/R_matrix.log','rt') as f:
                    for l in f:
                        print(l.strip())
        
        if sparse: R=csc_matrix(R)

        return R,ext




    def M(self,Te,ne,Ti=None,ni=None,E=0.1,sparse=False,write=True):
        from scipy.sparse import csc_matrix
        # Returns the rate matrix at temperature Te and density ne
        # Allow for different ion and electron parameters: if not defined,set to equal
        if Ti is None:
            Ti=Te
        if ni is None:
            ni=ne

        M,ext=self.populate('M',Te,ne,Ti,ni,E)
        if write:
            self.write_matrix(M,ext,'M',Te,ne,Ti,ne,E)
            if self.verbose:
                with open('logs/M_matrix.log','rt') as f:
                    for l in f:
                        print(l.strip())
        
        if sparse: M=csc_matrix(M)
    
        return M,ext
 
    def Se(self,Te,ne,Ti=None,ni=None,E=0.1):
        # Returns the energy loss rate matrix at temperature Te and density ne
        # Allow for different ion and electron parameters: if not defined,set to equal
        if Ti is None:
            Ti=Te
        if ni is None:
            ni=ne


    def gl_crm(self,Te,ne,Ti=None,ni=None,E=0.1,Sext=True,n=None):
        ''' Returns the P-space matrices '''
        from numpy import matmul,diag,real
        from numpy.linalg import inv,eig

        # Assume Te=Ti unless otherwise specified
        if Ti is None:
            Ti=Te
        # Assume ni=ne unless otherwise specified
        if ni is None:
            ni=ne
        # Assume intial distribution is defined in input unless otherwise specified
        if n is None:
            n=self.n0
        
        mat,ext=self.M(Te,ne,Ti,ni,E,write=False)

        # Create block matric from M
        MP=mat[:self.Np,:self.Np]
        MQ=mat[self.Np:,self.Np:]
        V=mat[self.Np:,:self.Np]
        H=mat[:self.Np,self.Np:]
        # Calculate Meff
        Meff=(MP-matmul(matmul(H,inv(MQ)),V))

        # Diagonalize M
        eigs,T=eig(mat)
        # Order the eigenvalues and vectors in increasing magnitude
        eigind=eigs.argsort()[::-1] 
        eigs=eigs[eigind]
        T=T[:,eigind]
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

        return Meff,GPp,nP0p


    def dndt(self,t,n,mat,ext):
        ''' Returns the time-derivative of the density for the CRM '''
        from numpy import matmul
        return matmul(mat,n)+ext


    def gl_nt(self,Te,ne,t,Ti=None,ni=None,E=0.1,n=None,Sext=True):
        ''' Solves the Greenland NpxNp problem '''
        if n is None:
            n=self.n0
        
        Meff,GPp,nP0p=self.gl_crm(Te,ne,Ti,ni,E,Sext,n)        

        return solve_ivp(lambda x,y: self.dndt(x,y,Meff,GPp),(0,t),nP0p,method='LSODA')


    def full_nt(self,Te,ne,t,n=None,Sext=True):
        ''' Solves the full NxN problem '''
        from scipy.integrate import solve_ivp,LSODA 
        if n is None:
            n=self.n0

        mat,ext=self.M(Te,ne,write=False)
        ext=(Sext is True)*ext
        return solve_ivp(lambda x,y: self.dndt(x,y,mat,ext),(0,t),self.n0,method='LSODA')



    def totparticles(self,arr,V=1):
        ''' Calculates the total particles of a density array with the same shape as self.species
            
            Assumes the order to be the same as in species
        '''
        from numpy import zeros
        try:
            ret=zeros((arr.shape[1],))
            for j in range(arr.shape[1]):
                for i in range(len(arr)):
                    ret[j]+=(1+('H2' in self.species[i]))*arr[i,j]
        except:
            ret=0
            for i in range(len(arr)):
                ret+=(1+('H2' in self.species[i]))*arr[i]

        return ret
        

    def create_UE_rates(self,fname='uerates',E=0.1,Sext=True):
        ''' Script that createst UEDGE rates '''
        from numpy import zeros
        from os import mkdir
        ret=zeros((15,60,self.Np,self.Np))
        ext=zeros((15,60,self.Np))
        for i in range(15):
            print('{}/15'.format(i+1))
            for j in range(60):
                ret[i,j,:,:],ext[i,j,:],_=self.gl_crm(10**(-1.2 +j/10),10**(10+0.5*i),E=E,Sext=Sext)

        # Assume species order to be [H0,H2]
        try:
            mkdir('output')
        except:
            pass
        
        with open('output/{}.dat'.format(fname),'w') as f:
            # H0 depletion
            f.write(' H0 depl. Rate(jt,jn) (cm**3 s**-1)  Te(jt) = 10**(-1.2 + (jt-1)/10) jt= 1,60\n')
            for j in range(15):
                f.write(' jn =   {}'.format(j+1)+(j==0)*'; jt = 1 -> 60 going by rows   ne(jn) = 10**(10 + 0.5*(jn-1)) jn=1,15'+'\n')
                out=''
                for k in range(10):
                    for i in range(6):
                        out+='{:1.5E}'.format(ret[j,i*10+k,0,0]).rjust(13,' ')
                    out+='\n'
                f.write(out+'\n')
            # H0 source from H2
            f.write(' H2->H0 Rate(jt,jn) (cm**3 s**-1)\n')
            for j in range(15):
                f.write(' jn =   {}\n'.format(j+1))
                out=''
                for k in range(10):
                    for i in range(6):
                        out+='{:1.5E}'.format(ret[j,i*10+k,0,1]).rjust(13,' ')
                    out+='\n'
                f.write(out+'\n')
            # H2 depletion
            f.write(' H2 depl. Rate(jt,jn) (cm**3 s**-1)\n')
            for j in range(15):
                f.write(' jn =   {}\n'.format(j+1))
                out=''
                for k in range(10):
                    for i in range(6):
                        out+='{:1.5E}'.format(ret[j,i*10+k,1,1]).rjust(13,' ')
                    out+='\n'
                f.write(out+'\n')
            # H2 creation
            f.write(' H2 creation Rate(jt,jn) (cm**3 s**-1)\n')
            for j in range(15):
                f.write(' jn =   {}\n'.format(j+1))
                out=''
                for k in range(10):
                    for i in range(6):
                        out+='{:1.5E}'.format(ret[j,i*10+k,1,0]).rjust(13,' ')
                    out+='\n'
                f.write(out+'\n')
            # H0 external source
            f.write(' H0 external source Rate(jt,jn) (cm**-3 s**-1)\n')
            for j in range(15):
                f.write(' jn =   {}\n'.format(j+1))
                out=''
                for k in range(10):
                    for i in range(6):
                        out+='{:1.5E}'.format(ext[j,i*10+k,0]).rjust(13,' ')
                    out+='\n'
                f.write(out+'\n')
            # H2 external source
            f.write(' H2 external source Rate(jt,jn) (cm**-3 s**-1)\n')
            for j in range(15):
                f.write(' jn =   {}\n'.format(j+1))
                out=''
                for k in range(10):
                    for i in range(6):
                        out+='{:1.5E}'.format(ext[j,i*10+k,1]).rjust(13,' ')
                    out+='\n'
                f.write(out+'\n')
    

            
'''=====================================================================
========================================================================
==============================REACTION==================================
========================================================================
====================================================================='''
 
class REACTION:
    # REACTION class
    # Contains:
    #   PARAMETERS
    #   name:       reaction name/handle
    #   type:       type of fit
    #   
    #   -energy differences
    #   
    #   FUNCTIONS
    #   rate(T):    returns the rate at temperature T

    def __init__(self, name, database, reactants, fragments, coeffs=0,typ='x',Tarr=0):
        from numpy import ones,array
        self.name=name
        self.coeffs=coeffs
        self.database=database
        # TODO: Tidy up this mess
        self.raw_reactants=[r.strip() for r in reactants.copy()]
        self.raw_fragments=[f.strip() for f in fragments.copy()]
        self.reactants=[r.strip() for r in reactants]
        self.fragments=[f.strip() for f in fragments]
        self.r_mult=ones((len(reactants),))
        self.f_mult=ones((len(fragments),))
        self.coeffs=array(coeffs)
        self.type=typ
        self.Tarr=array(Tarr)
        

        # Separate coefficients from reactant species
        for i in range(len(self.reactants)):
            if '*' in self.reactants[i]:
                mul=self.reactants[i].split('*')
                self.f_mult[i]=float(reactants[i].split('*')[0])
                self.reactants[i]=reactants[i].split('*')[1].strip()

        # Separate coefficients from fragment species
        for i in range(len(self.fragments)):
            if '*' in self.fragments[i]:
                mul=self.fragments[i].split('*')
                self.f_mult[i]=float(fragments[i].split('*')[0])
                self.fragments[i]=fragments[i].split('*')[1].strip()


    
    def print_reaction(self):
        # TODO tidy up this mess
        ret='{}_{}: '.format(self.database,self.name)
        

        return('{}_{}: {} => {}'.format(
            self.database,self.name,
            str(self.raw_reactants).replace('[','').replace('\'','').replace(']','').replace(' ','').replace(',',' + '),
            str(self.raw_fragments).replace('[','').replace('\'','').replace(']','').replace(' ','').replace(',',' + ')
                                        )
                ) 
            
    def rate(self,Te,Ti,E=None,ne=None,omegaj=1):
        ''' Returns the rate of the reaction '''
        from numpy import log,exp,sqrt,pi,inf,array
        from scipy.interpolate import interp1d,interp2d

        # Find reactant species
        if 'e' in self.reactants:
            T=Te
        elif 'p' in self.reactants:
            T=Ti
        else:
            ''' Photoemission '''
            T=0

        T=max(T,0.5)    # Validity range of temperature data according to Janev

        if self.type=='RATE':
            ''' We have an AMJUEL/similar rate '''
            ret=0
            if len(self.coeffs.shape)==2:
                ''' T,E fit '''
                for i in range(9):
                    for j in range(9):
                        ret+=self.coeffs[i,j]*(log(T)**i)*(log(E)**j)

            elif len(self.coeffs.shape)==1:
                ''' T fit '''
                for i in range(9):
                    ret+=self.coeffs[i]*(log(T)**i)
                    
            else:
                print('Unknown fit')

            return exp(ret)

        elif self.type=='COEFFICIENT':
            ''' Coefficient '''
            return self.coeffs

        elif self.type=='SIGMA':
            from scipy.integrate import quad
            # TODO Extend to general species?
            mH2=2*1.6735575e-27
            me=9.10938356e-31
            ev=1.602e-19
            VH2=sqrt((2*E*ev)/mH2)
            Ve=sqrt((2*T*ev)/me)
            mr=(me*mH2)/(me+mH2)
            vth=sqrt((2*self.coeffs[0]*ev)/mr)
            
            def sigma(E,Eth,q0,A,Omega,W,gamma,nu):
                Psi=(nu!=0)*(1-W/E)**nu+(gamma!=0)*(1-(W/E)**gamma)
                return (E>=Eth)*q0*(A/W**2)*((W/Eth)**Omega)*Psi
            
            def R(x,T,coeffs):
                me=9.10938356e-31
                return x*sigma(x*T,*self.coeffs)*exp(-x)

            return (4/sqrt(pi))*sqrt((T*ev)/(2*me))*quad(R,0,inf,args=(T,self.coeffs))[0]


        elif self.type=='ADAS':
            # TODO: How to deal with extrapolation!
            
            #f=interp1d(self.Tarr,self.coeffs,kind='slinear')
            Tu=list(self.Tarr)+[200,500,1e3,5e3,1e4,5e4,64000]
            su=list(self.coeffs)+[self.coeffs[-1],self.coeffs[-1],self.coeffs[-1],self.coeffs[-1],self.coeffs[-1],self.coeffs[-1],self.coeffs[-1]]
            f=interp1d(Tu,su,kind='slinear')

            # TODO: figure out what is implied by the statistical weight omegaj - set =1 for now
            return 2.1716e-8*(1/omegaj)*sqrt(13.6048/T)*f(T) 
        
        elif self.type=='UE':
            (jt,jn)=self.coeffs.shape
            ret=0
            t,n=[],[]
            
            for jt in range(jt):
                t.append(10**(-1.2+jt/10))
            for jn in range(jn):
                n.append(10**(10+0.5*jn))
            ne=min(n[-1],max(n[0],ne))

            f=interp2d(n,t,self.coeffs)
            return f(ne,T)[0]
                    
        else:
            print('Unknown type "{}"'.format(self.type))
    



'''=====================================================================
========================================================================
====================COLLISIONAL-RADIATIVE MODEL=========================
========================================================================
====================================================================='''

# Begin main function
from matplotlib.pyplot import figure
from numpy import linspace,amax
from scipy.integrate import solve_ivp
from time import sleep
crm=CRM()

'''
crm.create_UE_rates()
'''
ue=CRM(fname='input_simple.dat')
tf=5e-5
Te,ne=30,1e13
#nt,nttot=crm.nt(t,Te,ne)
#ntg,ntgtot=crm.ntg(t,Te,ne)


nt=crm.full_nt(Te,ne,tf,n=None,Sext=True)
ngt=crm.gl_nt(Te,ne,tf,n=None,Sext=True)
nue=ue.full_nt(Te,ne,tf,n=None,Sext=True)

#while lsoda.status=='running':
#    sleep(1)
#    print(lsoda.t)
#print(lsoda.status)


color=['b','r','g','c','m','gold','orange','grey','brown','lime']
f=figure()
ax=f.add_subplot(111)
ax.set_xlim((0,tf))
ntgtot,nttot,nuetot=crm.totparticles(ngt.y),crm.totparticles(nt.y),ue.totparticles(nue.y)

for i in range(crm.Np):
    ax.plot(nt.t,nt.y[i,:],linewidth=1.5,color=color[i],label=crm.species[i]+' full')
    ax.plot(ngt.t,ngt.y[i,:],'--',linewidth=3,label=crm.species[i]+' gl',color=color[i])
    ax.plot(nue.t,nue.y[i,:],':',linewidth=3,label=crm.species[i]+' UE',color=color[i])
for i in range(crm.Np,len(crm.species)):
    ax.plot(nt.t,nt.y[i,:])
ax.plot(nt.t,nttot,'k-',linewidth=1.5,label='Total full')
ax.plot(ngt.t,ntgtot,'k--',linewidth=3,label='Total gl')
ax.plot(nue.t,nuetot,'k:',linewidth=3,label='Total UE')
ax.legend(ncol=3)
ax.set_ylim((0,1.02*max([amax(ntgtot),amax(nttot),amax(nuetot)])))
ax.axhline(crm.totparticles(crm.n0),color='k',linewidth=0.5) 




        




