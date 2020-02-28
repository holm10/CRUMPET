# CRUM collisional-radiative model class reaction.py
# Separated from CRUM.py by holm10
# Changelog
# 200205 - Separated from CRUM.py #holm10

class CRM:
    def __init__(self,species,reactions,settings,path='.',recrad=None,ionizrad=None):
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
        self.reactions=reactions
        self.verbose=settings[0]
        self.Np=settings[1]
        self.n0=settings[2]
        self.path=path
        self.ionizrad=ionizrad
        self.recrad=recrad

        # Ensure that there is a logs directory under the run path
        try:
            mkdir('{}/logs'.format(self.path))
        except:
            pass


        # Write a log of the CRM setup path/logs
        with open('logs/setup.log','w') as f:
            f.write('CRUM run in {} on {}\n'.format(getcwd(),str(datetime.now())[:-7]))
            f.write('Defined species:\n')
            for i in self.species:
                f.write('    {}\n'.format(i))

            f.write('Defined reactions:\n')
            for r in self.reactions:
                f.write('{}\n'.format(r.print_reaction()))
        # Output to stdout if run verbosely
        if self.verbose:
            with open('logs/setup.log','rt') as f:
                for l in f:
                    print(l.strip())
        # Do the same for a Diagnostic rate matrix displaying reaction correlations
        self.DIAGNOSTIC()

        
    def getS(self,val,reactants,rad,Te,Ti,Tm,E,ne,bg):
        if isinstance(val,str):  
            temp=val.replace('erl1','0').replace('erl2','0')
            temp=temp.replace('Te',str(Te))
            temp=temp.replace('Ti',str(Ti))
            temp=temp.replace('Ta',str(Ti))         
            temp=temp.replace('Tm',str((Tm is not False)*Tm+(Tm is False)*E))

            
            ext=val.replace('erl1',str(rad*self.ionizrad.rate(Te,Ti,E,ne)))
            ext=ext.replace('erl2',str(rad*self.recrad.rate(Te,Ti,E,ne)))

            S=(bg in reactants)*eval(temp)

        else:
            S=(bg in reactants)*val  
        
        return(S)

        


    def populate(self,mode,Te,ne,Ti=None,ni=None,E=0,rad=True,Sind=None,Tm=False):
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
        from numpy import zeros
        
        if mode=='diagnostic':
            # Setup a 2D diagnostic list for the matrix and a list for the external source
            ext_source=[]
            ret=[]
            for i in range(len(self.species)):
                ext_source.append([])
                ret.append([])
                for j in range(len(self.species)):
                    ret[i].append([])
        else:
            # Setup a matrix and vector
            ret=zeros((len(self.species),len(self.species)))
            ext_source=zeros((len(self.species),))
 
        for i in range(len(self.species)):
            ''' Walk through each row (species)'''

            for r in self.reactions:
                ''' Sort the species of each reaction into the appropriate column '''
        
                ''' Get the energy loss term '''
                if mode=='S':
                    if Sind=='el':   # Calc electron loss
                        S=self.getS(r.S_r,r.reactants,rad,Te,Ti,Tm,E,ne,'e') 
                    elif Sind=='eg': # Electron interaction radiation
                        S=self.getS(r.S_g,r.reactants,rad,Te,Ti,Tm,E,ne,'e') 
                    elif Sind=='eV': # Electron interaction potential
                        S=self.getS(r.S_V,r.reactants,rad,Te,Ti,Tm,E,ne,'e') 
                    elif Sind=='pe': # Proton electron energy transfer
                        S=self.getS(r.S_e,r.reactants,rad,Te,Ti,Tm,E,ne,'p') 
                    elif Sind=='pg': # Proton interaction radiation
                        S=self.getS(r.S_g,r.reactants,rad,Te,Ti,Tm,E,ne,'p') 
                    elif Sind=='pV': # Proton interaction potential
                        S=self.getS(r.S_V,r.reactants,rad,Te,Ti,Tm,E,ne,'p') 
                
                


                
                # TODO: what if three-particle reaction?
                bg=('e' in r.reactants)*ne+('p' in r.reactants)*ni # Specify density for reactions
                if mode!='diagnostic':
                    bgm=(('e' in r.reactants)*ne)*(('p' in r.reactants)*ni) # Specify density for external source
                    bg=max(bg,1)    # Assure that auto-processes are considered
                j=None # Set flag to identify external sources

                # Loop through each reaction defined in the CRM
                for rea in range(len(r.reactants)):
                    # Find the column into which the fragments goes: if background mark external
                    try:
                        j=self.species.index(r.reactants[rea])  # Get the product species index of the correct column
                    except:
                        continue
            
                    multiplier=r.r_mult[rea] # Get the reaction multiplier for the density

                    if self.species[i]==r.reactants[rea]:   # If the species (row index) is a reactant, r is a depletion process 
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

                        elif mode=='S':
                            ''' Energy source matrix '''
                            # TODO: no self-sink?
                            ret[i,i]+=multiplier*r.rate(Te,Ti,E,ne)*bg*0 # Calculate the rate and store appropriately

                for frag in range(len(r.fragments)):    # Loop through the reaction fragments
                    ''' SOURCE '''
                    multiplier=r.f_mult[frag]   # Fragment multiplier

                    # Do nothing if background fragment
                    if r.fragments[frag] not in self.species: 
                        continue
                    
                    # If fragment enters row, add in appropriate column as defined by reactant
                    elif self.species.index(r.fragments[frag])==i:
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

                                elif mode=='S':
                                    ''' Energy source matrix '''
                                    ext_source[i]+=multiplier*r.rate(Te,Ti,E,ne)*bgm*S

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

                                elif mode=='S':
                                    ''' Energy loss matrix '''
                                    ret[i,j]+=multiplier*r.rate(Te,Ti,E,ne)*bg*S


        return ret,ext_source

        
    def write_matrix(self,mat,ext,char,te,ne,ti,ni,E):
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

        with open('logs/{}_matrix.log'.format(char),'w') as f:  # Create a file to write to according to char
            
            if char=='R': # Rate coefficient  matrix is being written
                f.write('Diagnostic rate coefficient (density-independend) matrix for CRUM run in {} on {}\n'.format(getcwd(),str(datetime.now())[:-7]))
                f.write('Te={} eV, Ti={} eV, E={} eV\n'.format(te,ti,E))
            elif char=='M': # Rate matrix is being written
                f.write('Diagnostic rate (density-dependend) matrix for CRUM run in {} on {}\n'.format(getcwd(),str(datetime.now())[:-7]))
                f.write('Te={} eV, Ti={} eV, ne={} 1/cm**3, ni={} 1/cm**3, E={} eV\n'.format(te,ti,ne,ni,E))
            elif char=='S': # Rate matrix is being written
                f.write('Diagnostic energy loss (density-dependend) matrix for CRUM run in {} on {}\n'.format(getcwd(),str(datetime.now())[:-7]))
                f.write('Te={} eV, Ti={} eV, ne={} 1/cm**3, ni={} 1/cm**3, E={} eV\n'.format(te,ti,ne,ni,E))
            
            # Create header line
            out='{}-MAT|'.format(char.upper()).rjust(10) 
            for s in self.species:
                out+=s.rjust(10,' ')
            f.write(out+'\n'+'_'*(1+len(self.species))*10+'\n')

            # Loop through the matrix
            for l in range(len(mat)):
                out=self.species[l]+'|' # Create a tabulated file with the row species displayed
                out=out.rjust(10,' ')
                # Add each element to the line with 10 characters reserved
                for e in mat[l,:]:
                    if e==0:
                        out+=' '*9+'_'
                    else:
                        out+='{:1.1E}'.format(e).rjust(10,' ')
                f.write(out+'\n')
            # Write the external source to the bottom of the output
            f.write('_'*(1+len(self.species))*10+'\n')
            out='S_ext|'.rjust(10,' ')
            for s in ext:
                out+='{:1.1E}'.format(s).rjust(10,' ')
            f.write(out)



    def DIAGNOSTIC(self):
        ''' Returns a diagnostic matrix containing lists of reactions accounted for in each element '''
        from os import getcwd
        from datetime import datetime
 
        dia,ext=self.populate('diagnostic',0,'*[ne]',0,'*[ni]') # Get the 2D diagnostic list and external source list

        # Write the diagnostic matrix to the logs
        with open('logs/reaction_matrix.log','w') as f:
            # Write the header line
            f.write('Diagnostic reaction matrix for CRUM run in {} on {}\n'.format(getcwd(),str(datetime.now())[:-7]))
            # Loop through each species
            for i in range(len(dia)):
                # Write the species header
                f.write('\n\n\n======== {} ========\n'.format(self.species[i]))
                # Start with the sinks
                f.write('{} DEPLETION:\n                 {}\n'.format(self.species[i],dia[i][i]))
                # Then do the sources
                for j in range(len(dia[i])):
                    if len(dia[i][j])>0:    # Don't write anything empty
                        if i!=j:    # Don't duplicate depletion
                            f.write('From {}:\n                 {}\n'.format(self.species[j],dia[i][j]))
                if len(ext[i])>0:   # Write external source last, if applicable
                    f.write('{} EXTERNAL SOURCE:\n                 {}\n'.format(self.species[i],ext[i]))
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
 
    def S(self,Te,ne,Ti=None,ni=None,E=0.1,rad=True,Tm=False,write=False):
        # TODO
        if Ti is None: Ti=Te # Check for Ti, set if necessary
        if ni is None: ni=ne # Check for ni, set if necessary

        Sel,extel=self.populate('S',Te,ne,Ti,ni,E,Sind='el',Tm=Tm)
        Seg,exteg=self.populate('S',Te,ne,Ti,ni,E,Sind='eg',Tm=Tm)
        SeV,extev=self.populate('S',Te,ne,Ti,ni,E,Sind='eV',Tm=Tm)
        Spe,extpe=self.populate('S',Te,ne,Ti,ni,E,Sind='pe',Tm=Tm)
        Spg,extpg=self.populate('S',Te,ne,Ti,ni,E,Sind='pg',Tm=Tm)
        SpV,extpv=self.populate('S',Te,ne,Ti,ni,E,Sind='pV',Tm=Tm)
        
        if write:
            self.write_matrix(Sel,extel,'Sel',Te,ne,Ti,ni,E)
            self.write_matrix(Seg,extel,'Seg',Te,ne,Ti,ni,E)
            self.write_matrix(SeV,extel,'SeV',Te,ne,Ti,ni,E)
            self.write_matrix(Spe,extel,'Spe',Te,ne,Ti,ni,E)
            self.write_matrix(Spg,extel,'Spg',Te,ne,Ti,ni,E)
            self.write_matrix(SpV,extel,'Spg',Te,ne,Ti,ni,E)
        
        return [Sel,Seg,SeV,Spe,Spg,SpV],[extel,exteg,extev,extpe,extpg,extpv]
    

    def S_gl(self,Te,ne,Ti=None,ni=None,E=0.1,rad=True,Sext=True,n=False,Tm=False,write=False):
        from numpy import matmul,concatenate
        from numpy.linalg import inv


        mat,ext=self.M(Te,ne,Ti,ni,E,write=write) # Get the full rate matrix
        
        U,Uext=self.S(Te,ne,Ti,ni,E,rad,Tm)
        

        if n is None: n=self.n0 # Use input n0 as default unless explict n0 requested
        
        # Create block matrix from M
        MP=mat[:self.Np,:self.Np]
        MQ=mat[self.Np:,self.Np:]
        V=mat[self.Np:,:self.Np]
        H=mat[:self.Np,self.Np:]
        nP0p=n[:self.Np]-matmul(matmul(H,inv(MQ)),n[self.Np:])

        [Sel,Seg,SeV,Spe,Spg,SpV],[extel,exteg,extev,extpe,extpg,extpv]=self.S(Te,ne,Ti,ni,E,rad,Tm)
        

        ret=[]
        for S in [  [Sel+Spe, extel+extpe],                     # e- energy source
                    [Spe+Spg+SpV-Sel, extpe+extpg+extpv-extel], # Ion/atom energy source
                    [Seg+Spg, exteg+extpg],                     # Radiation source
                    [SeV+SpV, extev+extpv]                      # Potential source
                 ]:
        


            UP=S[0][:self.Np,:self.Np]
            UQ=S[0][self.Np:,self.Np:]
            UV=S[0][self.Np:,:self.Np]
            UH=S[0][:self.Np,self.Np:]
        
            P=UP-matmul(UH,matmul(inv(MQ),V))
            Q=UV-matmul(UQ,matmul(inv(MQ),V))


            ret.append([concatenate((P,Q),axis=0),S[1],nP0p])
        
        return ret
        
    def gl_Et(self,Te,ne,t,Ti=None,ni=None,E=0.1,n=None,Tm=False,rad=True,Sext=True,write=False):
        from scipy.integrate import solve_ivp 
        from numpy import pad,array

        def dEdt(t,n,m1,m2,ext): return [m1*n[0]+m2*n[1]+ext,0]

        if n is None: n=self.n0 # Use input n0 as default unless explict n0 requested
        
        [Sel,Sia,Sg,SV]=self.S_gl(Te,ne,Ti,ni,E,rad,Sext,n,Tm,write) # Set up Greenland model 

        Selt,Siat,Sgt,SVt=[],[],[],[]
        for i in range(len(Sel[0])):
            Selt.append(solve_ivp(lambda x,y: dEdt(x,y,Sel[0][i,0],Sel[0][i,1],Sel[1][i]),(0,t),Sel[2]))#,method='LSODA'))
            Siat.append(solve_ivp(lambda x,y: dEdt(x,y,Sia[0][i,0],Sia[0][i,1],Sia[1][i]),(0,t),Sia[2]))#,method='LSODA'))
            Sgt.append(solve_ivp(lambda x,y: dEdt(x,y,Sg[0][i,0],Sg[0][i,1],Sg[1][i]),(0,t),Sg[2]))#,method='LSODA'))
            SVt.append(solve_ivp(lambda x,y: dEdt(x,y,SV[0][i,0],SV[0][i,1],SV[1][i]),(0,t),SV[2]))#,method='LSODA'))

        ret=[]
        for S in [Selt,Siat,Sgt,SVt]:
            s=[]
            for i in range(len(Sel[0])):
                s.append([array(S[i].t),array(S[i].y[0])])
            ret.append(S)

        return ret
        
       # return  [   solve_ivp(lambda x,y: self.dEdt(x,y,Sel[0],Sel[1]),(0,t),Sel[2],'LSODA'),
#                    solve_ivp(lambda x,y: self.dndt(x,y,Sia[0],Sia[1]),(0,t),Sia[2],method='LSODA'),
#                    solve_ivp(lambda x,y: self.dndt(x,y,Sg[0],Sg[1]),(0,t),Sg[2],method='LSODA'),
#                    solve_ivp(lambda x,y: self.dndt(x,y,SV[0],SV[1]),(0,t),SV[2],method='LSODA')
       #         ]

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

        if n is None: n=self.n0 # Use input n0 as default unless explict n0 requested
        

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
        eigind=abs(eigs).argsort()[::1] 
        eigs=eigs[eigind]
        T=T[:,eigind]
        # Only use real parts TODO is this OK??
        eigs=real(eigs)
        D=diag(eigs)
        if matrices is True:
            return mat,T,eigs

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

    def dEdt(self,t,n,mat,ext):
        ''' Returns the time-derivative of the density for the CRM i
            dndt(t,n,mat,ext)
        
            t   -   Time of evaluation [s]
            n   -   Vector of initial density distribution [cm**-3]
            mat -   Rate matrix
            ext -   External source vector
        '''
        from numpy import matmul,zeros
        ret=zeros((len(mat),))

        print((matmul(mat,n)+ext).shape)
        ret= matmul(mat,n)+ext
        return ret



    def dndt(self,t,n,mat,ext):
        ''' Returns the time-derivative of the density for the CRM i
            dndt(t,n,mat,ext)
        
            t   -   Time of evaluation [s]
            n   -   Vector of initial density distribution [cm**-3]
            mat -   Rate matrix
            ext -   External source vector
        '''
        from numpy import matmul

        return matmul(mat,n)+ext


    def gl_nt(self,Te,ne,t,Ti=None,ni=None,E=0.1,n=None,Sext=True):
        ''' Solves the Greenland NpxNp problem 
            gl_nt(Te,ne,t,*keys)
            
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

        if n is None: n=self.n0 # Use input n0 as default unless explict n0 requested
        
        Meff,GPp,nP0p=self.gl_crm(*self.M(Te,ne,Ti,ni,E,write=False),Sext,n) # Set up Greenland model 
        
        # Solve and return
        return solve_ivp(lambda x,y: self.dndt(x,y,Meff,GPp),(0,t),nP0p,method='LSODA')

    
    def generate_CRM(self,Te,ne,kappa,Ti=None,ni=None,E=0.1,n=None,Sext=True,epsilon=1):
        ''' Generates the optimal CRMs per Greenland 2001 '''
        from numpy import array,where,matmul
        from numpy.linalg import inv
        M,T,D=self.gl_crm(*self.M(Te,ne,Ti,ni,E,write=False),Sext=Sext,n=n,matrices=True) # Get matrices and eigenvalues/-vectors
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
        sorted_species=[self.species[i] for i in Zind]
        rej=[]
        for i in range(len(Z)):
            if Np[i]==True:
                if i+1==len(Z): continue # Full model, uninteresting
                TQ,delta=T[i+1:,i+1:],T[i+1:,:i+1]
                norm=matmul(inv(TQ),delta)
                maxnorm=max([sum(abs(norm[:,i])) for i in range(norm.shape[1])])
                if maxnorm<epsilon:
                    print('Possible CRM NP={} model with TQinv*delta={:.2E}, tauP={:.2E} s, tauQ={:.2E} s, and  P-space {} found.'.format(i+1,maxnorm,1/abs(min(sorted_D[:i+1])),1/abs(max(sorted_D[i+1:])),sorted_species[:i+1])) 
                else:
                    rej.append(maxnorm)
        if len(rej)>0:
            print('Minimum rejected maxnorm={:.2E}'.format(min(rej)))



    def full_nt(self,Te,ne,t,Ti=None,ni=None,E=0.1,n=None,Sext=True):
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

        if n is None: n=self.n0 # Use input n0 as default unless explict n0 requested

        mat,ext=self.M(Te,ne,Ti,ni,E,write=False) # Get the full rate matrix
        ext=(Sext is True)*ext  # Set source strength
        # Solve and return
        return solve_ivp(lambda x,y: self.dndt(x,y,mat,ext),(0,t),n,method='LSODA')



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
                    ret[j]+=(1+('H2' in self.species[i]))*arr[i,j]
        except: # If vector
            ret=0
            for i in range(len(arr)): # Loop through each species in arr
                ret+=(1+('H2' in self.species[i]))*arr[i]

        return ret
        

