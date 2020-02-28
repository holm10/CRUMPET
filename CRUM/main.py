# CRUM input parser class readinput.py
# Separated from CRUM.py by holm10
# Changelog
# 200205 - Separated from CRUM.py #holm10


class CRUM:
    def __init__(self,fname='input/CRUM.dat',path='.',vmax=14,nmax=8,verbose=False,NP=2):
        from CRUM.ratedata import RATE_DATA
        from CRUM.reactions import REACTION
        from CRUM.crm import CRM
        from numpy import zeros
        from os import getcwd
        ''' Generates the CRUM collisional-radiative model
        CRUM(fname='input/CRUM.dat',path='.')

        Looks for a CRUM input file in the location specified by path/fname
        
        Optional parameters
        vmax (14)   -   Number of excited vibrational molecular states considered
        nmax (8)    -   Number of atomic electronic levels considered
        verbose (F) -   Display more output
        Np (2)      -   P-space size, chosen as the first Np entries of the 'SPECIES' card
        '''
        self.path=path # Path to CRUM case




        def file2list(path,fname):
            ''' Returns a list of lines based on fname with comment and empty lines removed as well as the card and subcard locations in the list'''

            ret=[]
            # Parse out comments and empty lines and store data to memory
            with open('{}/{}'.format(path,fname)) as f:  # Open the file
                for l in f:
                    if len(l.split())==0: # Ignores empty lines
                        continue
                    elif l.strip()[0]=='#': # Lines starting with hashes are considered comments and ignored
                        continue
                    elif '#' in l:
                        ret.append(l.split('#')[0]) # Only store data appearing before a hash

                    else:   # Data lines are appended to list
                        ret.append(l)
        
            # Locate the cards in the list (separated by '**')
            cards=[i for i,x in enumerate(ret) if x[0:2]=='**']
            cards.append(len(ret))  # Append the end of the list for searching between entries
            # Locate the subcards in the list (separated by '*')
            subcards=[i for i,x in enumerate(ret) if x[0]=='*'] 
            subcards.append(len(ret))  # Append the end of the list for searching between entries

            return ret,cards,subcards # Return the list of lines, the lists of lines containing the starts of the cards and subcards




        def XY2num(string,X=None,Y=None):
            ''' Replaces 'X' or 'Y' with X and Y, respectively - or both, depending on which is defined '''

            if (X is None) and (Y is None): # No changes requested, return input
                return string
            elif Y is None:   # If no Y-values, only replace X
                    return string.replace('X',str(X))
            elif X is None: # If no X-values, only replace y
                    return string.replace('Y',str(Y))
            else:   # Both X and Y are replaced
                return string.replace('X',str(X)).replace('Y',str(Y))




        def rf(string):
            ''' Split reaction into reactants and fragments '''

            return string.split(' > ')


        def getS(energies,ene,X=None,Y=None):
            sauces=[k[:3] for k in energies]
            if 'S_e' not in sauces: S_e=None

            dia=[None,None,None,None]
            ener_s=[None,None,None,None]

            ret=[0,0,0,0]

            for i in energies:
                st=XY2num(i,X,Y)
            if st[:3]=='S_r': ener_s[0]=st
            elif st[:3]=='S_g': ener_s[1]=st
            elif st[:3]=='S_V': ener_s[2]=st
            elif st[:3]=='S_e': ener_s[3]=st


            for i in energies:
                st=XY2num(i,X,Y)

                if i[:3] in ['S_r','S_g','S_V','S_e']:
                    

                    for v in list(ene.keys())[::-1]:
                        st=st.replace(str(v),str(ene[v]))

                        
                    
                    if i[:3]=='S_r': 
                        dia[0]=st
                        try:
                            exec(st,globals())
                            ret[0]=S_r
                        except:
                            ret[0]=st
                    elif i[:3]=='S_g':
                        dia[1]=st
                        try:
                            exec(st,globals())
                            ret[1]=S_g
                        except:
                            ret[1]=st
                    elif i[:3]=='S_V':
                        dia[2]=st
                        try:
                            exec(st,globals())
                            ret[2]=S_V
                        except:
                            ret[2]=st
                    elif i[:3]=='S_e':
                        dia[3]=st
                        try:
                            exec(st,globals())
                            ret[3]=S_e
                        except:
                            ret[3]=st

             
            #print(energies)
            #print(ener_s)
            #print(dia)
            #print(ret)

            return ret


        # Read the input file into a list
        lines,cards,subcards=file2list(path,fname)



                
        ''' LOOP THROUGH THE DEFINED CARDS '''
        for i in range(len(cards)-1):

            ''' Store species from species card to temporary list '''
            if lines[cards[i]].split()[1].upper()=='SPECIES':
                self.species=[]
                for j in range(cards[i]+1,cards[i+1]): # Look between this card and the next
                    self.species.append(lines[j].split()[0].strip()) # Add species to class object
                # Set up array for initial densities
                n0=zeros((len(self.species),))


                '''  Store the reactions requested to temporary reaction list '''
            elif lines[cards[i]].split()[1].upper()=='REACTIONS':
                reactionlist=[]
                # TODO: Correct once energies enter the input deck
#                for j in range(cards[i]+1,cards[i+1],2): # Look between this card and the next: two lines per reactions
                for k in range(subcards.index(cards[i])+1,subcards.index(cards[i+1])):

                    rname=lines[subcards[k]].split()[1]   # Name is the second entry on the first line: subcard first entry, discarded

                    if rname.upper()=='CUSTOM': # If custom reactions deck, store the filepath, which is the only entry on subcard
                        reactionlist.append(    [   rname, lines[subcards[k]+1].strip()   ]   ) # Store name and path

                    else: # If regular reaction, create two lists with reactants and fragments and unpack
                        energies=[]
                        for l in range(subcards[k]+2,subcards[k+1]): # Loop through the lines following the declaration looking for energy sinks/sources
                            energies.append(lines[l].strip())

                        reactionlist.append(    [   rname, *[x.split(' + ') for x in rf(lines[subcards[k]+1])], energies  ]   ) # Store name, reactants and framgents
                        


                ''' Store the locations of the standard rate data files relative to CWD'''
            elif lines[cards[i]].split()[1].upper()=='RATES':
                for j in range(cards[i]+1,cards[i+1]): # Look between this card and the next
                    # Extract the database type and path
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
                # Create RATE_DATA class object based on the above files paths
                self.ratedata=RATE_DATA(amjuel,hydhel,h2vibr,ADAS,UE,getcwd())



                ''' Setup the CRM '''
            elif lines[cards[i]].split()[1].upper()=='SETTINGS':
                for j in range(cards[i]+1,cards[i+1]): # Look between this card and the next
                    # Create setting based on card type
                    try:
                        [c,setting,value]=lines[j].split()
                    except:
                        [c,setting]=lines[j].split()

                    setting=setting.upper() # Uppercase for comparative purposes
                    if c!='*':  # Ensure that we are reading subcards
                        continue
                    elif setting=='VMAX':
                        vmax=int(value)
                    elif setting=='NMAX':
                        nmax=int(value)
                    elif setting=='VERBOSE':
                        verbose=bool(int(value))
                    elif setting=='NP':
                        self.Np=int(value)
                    elif setting=='N0':
                        for k in range(j+1,subcards[subcards.index(j)+1]): # Loop through all defined initial densities
                            # Store the density in the location of n0 that corresponds to the species requested
                            n0[self.species.index(lines[k].split()[0])]=float(lines[k].split()[1]) 
                    else:
                        print('Unrecognized setting {}. Ignoring.'.format(setting))
                        continue


                ''' Store energies '''
            elif lines[cards[i]].split()[1].upper()=='ENERGIES':
                elst,_,_=file2list(path,lines[cards[i]+1].strip())
                ene=dict() # Dictionary for energy handles
                for el in elst: # Loop through each line
                    ene[el.split('=')[0].strip()]=float(el.split('=')[1].strip()) 

                # Make energy handles available within function
                locals().update(ene) 


                ''' Unknown card, abort '''
            else:
                print('Unknown card "{}": aborting!'.format(lines[cards[i]].split()[1]))
                return


         

        ''' END CARDS LOOP '''

        # Create list of reaction objects based on the previously read reaction list
        reactions=[]

        ''' LOOP OVER ALL THE DEFINED REACTIONS '''
        for r in reactionlist:



            # Split the definitions into names and databases
            database=r[0].split('_')[0]
            try:
                ID=r[0].split('_')[1]
            except:
                pass

            ''' Vibrational transitions in molecules '''
            if 'XvY' in ID:
                for x in range(vmax+1): # Loop over the requested excited vibrational states
                    for y in [-1,1]: # Assume only ladder-like vibrational transitions (dv=+/-1)
                        if x+y in range(vmax+1): # Make sure we do not excite beyond vmax
                            vID=XY2num(ID,x,x+y) # Substitute initial and final state into temporary string
                            
                    

                            # Create reaction object
                            reactions.append(REACTION(              vID, # Reaction name
                                                                    database, # Reaction database
                                                                    [XY2num(i,x,x+y) for i in r[1]], # Reactant with vib. level numbers
                                                                    [XY2num(i,x,x+y) for i in r[2]], # Product with vib. level numbers
                                                                    self.ratedata.get_coeff(database,vID), # Coefficients of the reaction
                                                                    'RATE', # Type of reaction
                                                                    getS(r[-1],ene,x,x+y)
                                                                    
                                                        )
                                            )
            
                ''' Vibrationally-dependent reaction '''
            elif 'Xl' in ID:
                for x in range(vmax+1): # Generate reactions for each vibrational level
                    vID=XY2num(ID,x,x) # Substitute initial and final state into temporary string
                    reactions.append(REACTION(              vID, # Reaction name
                                                            database, # Reaction database
                                                            [XY2num(i,x) for i in r[1]], # Reactant with vib. level numbers
                                                            [XY2num(i,x) for i in r[2]], # Product with vib. level numbers
                                                            self.ratedata.get_coeff(database,vID), # Coefficients of the reaction
                                                            'RATE', # Type of reaction
                                                            getS(r[-1],ene,x)

                                                )
                                    )



                ''' Custom reaction-rate '''
            elif database.upper()=='CUSTOM':
                # Open the custom data file and read lines into data and subcard locations into list
                data,_,subcards=file2list(path,r[1].strip())
                database_use=data[0].strip()    # Use the first line as database

                # Loop through the cards
                for i in range(len(subcards)-1):
                    subc=subcards[i] # Helper index
                    fit=data[subc].split()[1].upper()   # What kind of fit is being used - first subcard entry
                    name=data[subc].split()[2]  # What is the name of the reaction - second subcard entry
                    reactants,fragments=rf(data[subc+1])    # Get the reactants and fragments from the second line

                    

                    # Identify the type of reaction we are working with
                    ''' v-dependent rate '''
                    if fit=='RATE':
                        eng=dict()
                        if 'X' in name:
                            # TODO: what if not all vib. states have data?
                            for j in range(vmax+1): # Read data for each vibrational state up to vmax
                                m=0
                                engs=[]
                                while data[m][0]!='v':
                                    if data[m][0:2]=='S_':
                                        engs.append(data[m].strip()) 
                                    m+=1
                                    

                                reactions.append(REACTION(              XY2num(name,j), # Name
                                                                        database_use,   # Database
                                                                        XY2num(reactants,j).split(' + '), # Reactants w/ v-level number
                                                                        XY2num(fragments,j).split(' + '), # Products w/ v-level number
                                                                        [float(x) for x in data[subc+m+j*2].split()], # Coefficients for v-level
                                                                        'RATE',  # Type of reaction
                                                                        getS(engs,ene,j)
                                                                    )
                                                     )
                        # TODO: what if non-v dependent rate?

                        ''' Transition coefficient '''
                    elif fit=='COEFFICIENT':
                        if 'X' in name: # v-dependent rate
                            # TODO: what if not all vib. states have data?
                            m=2
                            engs=[]
                            while data[subc+m][:2]=='S_':
                                engs.append(data[subc+m].strip()) 
                                m+=1

                            for j in range(vmax+1): # Read data for each vibrational state up to vmax
                                    



                                reactions.append(REACTION(          XY2num(name,j), # Name
                                                                    database_use,   # Database
                                                                    XY2num(reactants,j).split(' + '), # Reactants w/ v-level number
                                                                    XY2num(fragments,j).split(' + '), # Products w/ v-level number
                                                                    float(data[subc+m+j].split()[1]), # Transition coefficient for v-level
                                                                    'COEFFICIENT', # Type of reaction
                                                                    getS(engs,ene,j)
                                                                   )
                                                     )
                        else: # Other transition
                            m=2
                            engs=[]
                            while data[subc+m][:2]=='S_':
                                engs.append(data[subc+m].strip()) 
                                m+=1
                                    
                            reactions.append(REACTION(              name,           # Name 
                                                                    database_use,   # Database
                                                                    reactants.split(' + '), # Reactants
                                                                    fragments.split(' + '), # Fragments
                                                                    float(data[subc+m]),    # Transition coefficient
                                                                    'COEFFICIENT',   # Type of reaction
                                                                    getS(engs,ene)
                                                                )
                                                 )

                        ''' Cross-section '''
                    elif fit=='SIGMA':
                        m=2
                        engs=[]
                        while data[subc+m][:2]=='S_':
                            engs.append(data[subc+m].strip()) 
                            m+=1
                        # TODO: extend definitions of cross-section
                        # Presently assumes SAWADA-like cross-section definition
                        reactions.append(REACTION(              name,           # Name
                                                                database_use,   # Database
                                                                reactants.split(' + '), # Reactants
                                                                fragments.split(' + '), # Fragments
                                                                [float(x) for x in data[subc+m].split()], # SAWADA cross-section parameters
                                                                'SIGMA', # Type of reaction
                                                                getS(engs,ene)
                                                            )
                                             )
                    else:
                        print('Reaction type "{}" not recognized! Aborting.'.format(data[subcards[i]].split()[1]))
                        return
                    

                ''' ADAS data '''
            elif database.upper()=='ADAS':
                for x in range(1,nmax+1): # Read the data for each electronic state up to nmax


                    if ID.upper()=='EXCITATION': # Electron impact excitation
                        rn=range(x+1,nmax+1) # Excitation only possible between current state and nmax
                        fit='ADAS'          # ADAS-type fit
                        Tarr=self.ratedata.get_coeff(database,'T')  # Temperature array for interpolation


                    elif ID.upper()=='RELAXATION': # Radiative relaxation
                        rn=range(1,x)   # Relaxation only possible to lower states
                        fit='COEFFICIENT' # Coefficient-like decay
                        Tarr=None   # Not T-dependent

                    # Loop through each of the available final states
                    for y in rn:
                            reactions.append(REACTION(              '{}_{}-{}'.format(ID,x,y),  # Name w/ appended initial-final state
                                                                    database,                   # Database
                                                                    [XY2num(i,x,y) for i in r[1]],  # Reactants 
                                                                    [XY2num(i,x,y) for i in r[2]],  # Fragments
                                                                    self.ratedata.get_coeff(database,'{}-{}'.format(x,y)),  # Get transition coefficient
                                                                    fit,    # Reaction type
                                                                    getS(r[-1],ene,x,y),
                                                                    Tarr    # Temperature array for interpolation
                                                                )
                                                 )


                ''' UEDGE rates '''
            elif database.upper()=='UE':
                reactions.append(REACTION(              ID,         # Name
                                                        database,   # Database
                                                        r[1],       # Reactants
                                                        r[2],       # Fragments
                                                        self.ratedata.get_coeff(database,ID),   # Get coefficients
                                                        'UE',        # Reaction type
                                                        getS(r[-1],ene)
                                                   )
                                     )
                            
                    
                ''' EIRENE rates '''
            elif database.upper() in ['HYDHEL','AMJUEL','H2VIBR']:




                reactions.append(REACTION(              ID,         # Name
                                                        database,   # Database
                                                        r[1],       # Reactants
                                                        r[2],       # Fragments
                                                        self.ratedata.get_coeff(database,ID), # Get coefficients
                                                        'RATE',      # Reaction type
                                                        getS(r[-1],ene)
                                                    )
                                     )
            
            else:
                print('Database "{}" not recognized! Aborting.'.format(database))
                return
 

            ''' END LOOP OVER DEFINED REACTIONS '''


        # Define reactions for UEDGE raditation
        ionizrad=REACTION('IONIZRAD','UE','','',self.ratedata.get_coeff('UE','IONIZRAD'),'UE',[None,None,None,None])
        recrad=REACTION('RECRAD','UE','','',self.ratedata.get_coeff('UE','RECRAD'),'UE',[None,None,None,None])


        # Setup the crm
        self.crm=CRM(self.species,reactions,[verbose,self.Np,n0],self.path,recrad=recrad,ionizrad=ionizrad)
        

    def totpart(self,arr,V=1):
        ''' Calculates the total particles in the array arr
            totpart(arr)

            Optional parameters
            V (1)   -   Volume of the assumed domain in cm**3

            Uses the CRM function totparticles
        '''
        return self.crm.totparticles(arr,V)


    def n_gl(self,Te,ne,t,Ti=None,ni=None,E=0.1,n=None,Sext=True):
        ''' Calculates the Greenland (P-space) density evolution in a 1cm**3 box up to t
            n_gl(Te,ne,t)
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
            
            Uses the CRM function gl_nt
        '''
        return self.crm.gl_nt(Te,ne,t,Ti,ni,E,n,Sext)

    
    def n_full(self,Te,ne,t,Ti=None,ni=None,E=0.1,n=None,Sext=True):
        ''' Calculates the full CRM density evolution in a 1cm**3 box up to t
            n_gl(Te,ne,t)
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
        return self.crm.full_nt(Te,ne,t,Ti,ni,E,n,Sext)

    def create_UE_rates(self,fname='uerates',E=0.1,Sext=True,h0h2=['H(n=1)','H2(v=0)']):
        ''' Script that writes UEDGE rates to self.path/fname.dat
            create_UE_rates(fname='uerates')

            Optional parameters
            E (0.1)     -   target particle energy [eV]
            Sext (True) -   Include external source (from background plasma reactions into CRM species)
            h0h2 (['H(n=1)','H2(v=0)']) - List of the H0 and H2 at their electronic and vibrational ground state handles used in the input
            '''
        from numpy import zeros
        from os import mkdir

        # Ensure that the CRM format coincides with the assumed format
        if self.species[:2]!=h0h2:
            print('This model assumes the first species to be H0 and the second to be H2 in their electronic and vibrational ground states.\nThe first input species do not match {} and {}, respectively. Aborting.'.format(h0h2[0],h0h2[1]))
            return

        if self.Np!=2:
            print('Presently, UEDGE only considers two neutral hydrogenic species, H0 and H2. The P-space of the model is larger than thisand, thus, is incompatible with the UEDGE model. No rates can be written for Np!=2. Aborting.')
            return

        # Create deinsty and temperature points in log-log space as in existing UEDGGE rate files
        ret=zeros((15,60,self.Np,self.Np))
        ext=zeros((15,60,self.Np))
        # Calculate the Greenland (Np) space rates for the T-n space
        print('Creating UEDGE rate output')
        for i in range(15): # 15 density points in log-log space
            print('    Density step {}/15'.format(i+1))
            for j in range(60): # 60 temperature points in log-log space
                ret[i,j,:,:],ext[i,j,:],_=self.crm.gl_crm(10**(-1.2 +j/10),10**(10+0.5*i),E=E,Sext=Sext) # Store to matrix

        # Ensure output directory exists 
        try:
            mkdir('{}/output'.format(self.path))
        except:
            pass
        
        # Write UEDGE data
        print(' Writing UEDGE rates to {}'.format(fname))
        with open('output/{}.dat'.format(fname),'w') as f:
            # Headers and positions in Greenland matrices
            setups= [   [' H0 depl. Rate(jt,jn) (cm**3 s**-1)  Te(jt) = 10**(-1.2 + (jt-1)/10) jt= 1,60\n', ret[:,:,0,0]    ],
                        [' H2->H0 Rate(jt,jn) (cm**3 s**-1)\n', ret[:,:,0,1]                                                ],
                        [' H2 depl. Rate(jt,jn) (cm**3 s**-1)\n', ret[:,:,1,1]                                              ],
                        [' H2 creation Rate(jt,jn) (cm**3 s**-1)\n', ret[:,:,1,0]                                           ],
                        [' H0 external source Rate(jt,jn) (cm**-3 s**-1)\n', ext[:,:,0]                                     ],
                        [' H2 external source Rate(jt,jn) (cm**-3 s**-1)\n', ext[:,:,1]                                     ]
                    ]
      
            for l in range(len(setups)): # Loop through all species contributions
                print('    Writing rate {}/{}'.format(l+1,len(setups)))
                f.write(setups[l][0])
                for j in range(15): # Write each density block
                    out=''
                    f.write(' jn =   {}'.format(j+1)+(j==0)*(l==0)*'; jt = 1 -> 60 going by rows   ne(jn) = 10**(10 + 0.5*(jn-1)) jn=1,15'+'\n')
                    for k in range(60): # Store temperature data to string  
                        out+='{:1.5E}'.format(setups[l][1][j,k]).rjust(13,' ')
                        if (k+1)/6 in range(1,11):  # Split string into ten rows
                            out+='\n'
                    f.write(out+'\n')   # Write the data to file


    def plot_Et(self,t,Te,ne,E=0.1,Ti=None,ni=None,Sext=True,Np=True,Nq=False,N=False,fig=None,n0=None,ax=0,figsize=(10,10/1.618033),linestyle='-',ylim=None,linewidth=2,labelapp='',qlabel=False,savename=None,title=None,pretitle='',figtype='png',color=None,ncol=3,nuclei=True,Tm=False,rad=True,gl=False):
        ''' Plots the time-evolution of the full CRM or Greenland model on a plasma background 
            plot_nt(t,Te,ne,*keys)

            t   -   Maximum simulations time to be plotted [s]
            Te  -   Background electron plasma temperature [eV]
            ne  -   Background electron plasma density [cm**-3]
            
            Optional parameters
            E (0.1)         -   Target particle energy [eV]
            Ti (None)       -   Background ion plasma temperature, Ti=Te assumed if None [eV]
            ni (None)       -   Background ion plasma density, ni=ne assumed if None [cm**-3]
            Sext (True)     -   Switch whether to consider sources from the background plasma or not
            gl (False)      -   Switch whether to plot the Greenland CRM (P-space only). 
                                If false, full CRM is evaluated, and plotting is controlled by Np, Nq and N.
            Np (True)       -   Switch whether to plot the P-space evolution if the full model is evaluated
            Nq (False)      -   Switch whether to plot the Q-space evolution if the full model is evaluated
            N (False)       -   Switch whether to plot all species if the full model is evaluated
            fig (None)      -   Figure object on which to plot. Creates new figure if None or unviable
            ax (0)          -   Which axis of fig to plot on, if fig is viable for plotting
            figsize (10,16) -   Figure size, if created by script
            linestyle ('-') -   Plot line style for the requested spaces
            linewidth (2)   -   Line width of P-space and all lines if N-space; Q-space linewidth=0.5*linewidth
            color (None)    -   User-defined list of colors for P-space plotting
            title (None)    -   User-defined title, overrides auto-generated title
            pretitle ('')   -   Prefix for auto-generated title (temperature and density of BG plasma)
            labelapp ('')   -   String to append to labels in legend
            qlabel (False)  -   Switch for displaying the Q-space labels on the legend
            ncol (3)        -   Number of columns in the legend
            ylim (None)     -   Requested y-limits: defaults to axis max if None
            savename (None) -   Figure savename. Default location path/output/figs. Not saved if None
            figtype ('png') -   Figure type if saved
            nuclei (True)   -   The total number of nuclei is plotted if True (2*n_m), the number of particles if False
            
        
        '''
        from matplotlib.pyplot import figure
        from numpy import log10,sum,matmul,zeros
        from numpy.linalg import inv
        from os import mkdir
        
        # Get axis handle or create figure
        try:
            ax=fig.get_axes()[0]
        except:
            fig=figure(figsize=figsize)
        
        ev=1.602e-19
        # If no initial distribution is requested, use the one specified in the input
        if n0 is None:
            n0=self.crm.n0

        # Check what model to use, Greenland or full
        if gl is True: # Greenland
            nt=self.n_gl(Te,ne,t,Ti,ni,E,n0,Sext) # Solve ODE
        else: # Full model
            nt=self.n_full(Te,ne,t,Ti,ni,E,n0,Sext) # Solve ODE

        Sel_p,Sia_p,Sg_p,SV_p=zeros(nt.y.shape),zeros(nt.y.shape),zeros(nt.y.shape),zeros(nt.y.shape),

        [Sel,Sia,Sg,SV]=self.crm.S_gl(Te,ne,Ti,ni,E,rad,Sext,self.crm.n0,Tm) # Set up Greenland model 


        mat,ext=self.crm.M(Te,ne,Ti,ni,E,write=False) # Get the full rate matrix
        # Create block matrix from M
        MP=mat[:self.Np,:self.Np]
        MQ=mat[self.Np:,self.Np:]
        V=mat[self.Np:,:self.Np]
        H=mat[:self.Np,self.Np:]
        # Manually integrate
        for dt in range(len(nt.t)-1):
            # Calculate the new inital distribution
            nP0p=nt.y[:self.Np,dt]-matmul(matmul(H,inv(MQ)),nt.y[self.Np:,dt])

            # Now, for each time-step, calculate the 
            Sel_p[:,dt]=(matmul(Sel[0],Sel[2])+Sel[1])*(nt.t[dt+1]-nt.t[dt])
            Sia_p[:,dt]=(matmul(Sia[0],Sel[2])+Sia[1])*(nt.t[dt+1]-nt.t[dt])
            Sg_p[:,dt]=(matmul(Sg[0],Sel[2])+Sg[1])*(nt.t[dt+1]-nt.t[dt])
            SV_p[:,dt]=(matmul(SV[0],Sel[2])+SV[1])*(nt.t[dt+1]-nt.t[dt])
                   
        if color is None: # Define color sequence up to 10 unless specific sequence requested
            color=[ 'b', 'r', 'm', 'c', 'darkgreen', 'gold', 'brown' ,'lime', 'grey', 'orange' ]


        pl=[-Sel_p,Sia_p,Sg_p,SV_p,Sel_p+Sia_p+Sg_p+SV_p]
        title=['Electron E loss', 'Ion/atom E source','Radiation','Binding energy','Power balance']
        
        for j in range(5):
            if j!=4:
                ax=fig.add_subplot(3,2,j+1)
            else:
                ax=fig.add_subplot(3,1,3)

            # Plot P-space species if requested
            if Np:
               for i in range(self.Np):
                    ax.plot(nt.t[:-1]*1e3,pl[j][i,:-1]*ev,linewidth=linewidth,color=color[i],label=self.species[i]+labelapp,linestyle=linestyle)
            # Plot Q-space species if requested
            if Nq:
                for i in range(self.Np+1,len(self.species)):
                    if qlabel:
                        ax.plot(nt.t[:-1]*1e3,pl[j][i,:-1]*ev,linewidth=linewidth*(0.5)**Np,label=self.species[i]+labelapp,linestyle=linestyle)
                    else:
                        ax.plot(nt.t[:-1]*1e3,pl[j][i,:-1]*ev,linewidth=linewidth*(0.5)**Np,linestyle=linestyle)
            # Plot all species w/ automated colors if requested
            if N:
               for i in range(len(self.species)):
                    if qlabel:
                        ax.plot(nt.t[:-1]*1e3,pl[j][i,:-1]*ev,linewidth=linewidth,label=self.species[i]+labelapp,linestyle=linestyle)
                    else:
                        ax.plot(nt.t[:-1]*1e3,pl[j][i,:-1]*ev,linewidth=linewidth,linestyle=linestyle)
        
            ax.set_title(title[j])
            ax.set_xlabel('Time [ms]') # Show x-label
            ax.set_ylabel(r'Power [W]') # Show y-label
            #ax.plot(nt.t[:-1]*1e3,self.totpart(pl[j])[:-1]*ev,'k',linewidth=linewidth,label='Total'+labelapp,linestyle=linestyle)



        # Unless title is requested, show automated title
        if title is None: # Automated title
            ex=int(log10(ne))
            mu=ne/(10**ex)
            ax.set_title(pretitle+r' $\rm{{T_e}}$={} eV, $\rm{{n_e}}$={}$\times10^{{{}}}$ $\rm{{cm^{{-3}}}}$'.format(Te,mu,ex))
        else: # Requested title
            ax.set_title(title)

        ax.legend(ncol=ncol,loc='best') # Show legend
        ax.set_xlim(0,t*1e3) # Set X-lim in ms
        # If no ylim is requested set the ylim to the highest plotted line
        '''
        if ylim is None:
            # Total particles in this plot
            totmax=max(self.totpart(nt.y))
            # Compare to previous max
            if totmax>ax.get_ylim()[1]:
                # If higher, update
                ax.set_ylim((0,totmax))
            else: # Set bottom to zero
                ax.set_ylim(bottom=0)
        else: # Use requested ylim
            ax.set_ylim(ylim)
        '''
        
        # Save to casepath/output/figs/figname.type if figname set
        if savename is not None:
            # Make sure that the direcory exists
            try:
                mkdir('{}/output/figs'.format(self.path))
            except:
                pass
            # Save
            fig.savefig('{}/output/figs/{}.{}'.format(self.path,savename,figtype),dpi=300,edgecolor=None,format=figtype,bbox_inches='tight')


        fig.show() # Show fig
        return fig # Return figure object
            
            



    def plot_nt(self,t,Te,ne,E=0.1,Ti=None,ni=None,Sext=True,Np=True,Nq=False,N=False,fig=None,n0=None,ax=0,figsize=(10,10/1.618033),linestyle='-',gl=False,ylim=None,linewidth=2,labelapp='',qlabel=False,savename=None,title=None,pretitle='',figtype='png',color=None,ncol=3,nuclei=True):
        ''' Plots the time-evolution of the full CRM or Greenland model on a plasma background 
            plot_nt(t,Te,ne,*keys)

            t   -   Maximum simulations time to be plotted [s]
            Te  -   Background electron plasma temperature [eV]
            ne  -   Background electron plasma density [cm**-3]
            
            Optional parameters
            E (0.1)         -   Target particle energy [eV]
            Ti (None)       -   Background ion plasma temperature, Ti=Te assumed if None [eV]
            ni (None)       -   Background ion plasma density, ni=ne assumed if None [cm**-3]
            Sext (True)     -   Switch whether to consider sources from the background plasma or not
            gl (False)      -   Switch whether to plot the Greenland CRM (P-space only). 
                                If false, full CRM is evaluated, and plotting is controlled by Np, Nq and N.
            Np (True)       -   Switch whether to plot the P-space evolution if the full model is evaluated
            Nq (False)      -   Switch whether to plot the Q-space evolution if the full model is evaluated
            N (False)       -   Switch whether to plot all species if the full model is evaluated
            fig (None)      -   Figure object on which to plot. Creates new figure if None or unviable
            ax (0)          -   Which axis of fig to plot on, if fig is viable for plotting
            figsize (10,16) -   Figure size, if created by script
            linestyle ('-') -   Plot line style for the requested spaces
            linewidth (2)   -   Line width of P-space and all lines if N-space; Q-space linewidth=0.5*linewidth
            color (None)    -   User-defined list of colors for P-space plotting
            title (None)    -   User-defined title, overrides auto-generated title
            pretitle ('')   -   Prefix for auto-generated title (temperature and density of BG plasma)
            labelapp ('')   -   String to append to labels in legend
            qlabel (False)  -   Switch for displaying the Q-space labels on the legend
            ncol (3)        -   Number of columns in the legend
            ylim (None)     -   Requested y-limits: defaults to axis max if None
            savename (None) -   Figure savename. Default location path/output/figs. Not saved if None
            figtype ('png') -   Figure type if saved
            nuclei (True)   -   The total number of nuclei is plotted if True (2*n_m), the number of particles if False
            
        
        '''
        from matplotlib.pyplot import figure
        from numpy import log10,sum
        from os import mkdir
        
        # Get axis handle or create figure
        try:
            ax=fig.get_axes()[0]
        except:
            fig=figure(figsize=figsize)
            ax=fig.add_subplot(111)
        
        # If no initial distribution is requested, use the one specified in the input
        if n0 is None:
            n0=self.crm.n0

        # Check what model to use, Greenland or full
        if gl is True: # Greenland
            nt=self.n_gl(Te,ne,t,Ti,ni,E,n0,Sext) # Solve ODE
            Np,Nq,N=True,False,False # Np only option
        else: # Full model
            nt=self.n_full(Te,ne,t,Ti,ni,E,n0,Sext) # Solve ODE

        if color is None: # Define color sequence up to 10 unless specific sequence requested
            color=[ 'b', 'r', 'm', 'c', 'darkgreen', 'gold', 'brown' ,'lime', 'grey', 'orange' ]

        # Plot P-space species if requested
        if Np:
           for i in range(self.Np):
                ax.plot(nt.t*1e3,nt.y[i,:],linewidth=linewidth,color=color[i],label=self.species[i]+labelapp,linestyle=linestyle)
        # Plot Q-space species if requested
        if Nq:
            for i in range(self.Np+1,len(self.species)):
                if qlabel:
                    ax.plot(nt.t*1e3,nt.y[i,:],linewidth=linewidth*(0.5)**Np,label=self.species[i]+labelapp,linestyle=linestyle)
                else:
                    ax.plot(nt.t*1e3,nt.y[i,:],linewidth=linewidth*(0.5)**Np,linestyle=linestyle)
        # Plot all species w/ automated colors if requested
        if N:
           for i in range(len(self.species)):
                if qlabel:
                    ax.plot(nt.t*1e3,nt.y[i,:],linewidth=linewidth,label=self.species[i]+labelapp,linestyle=linestyle)
                else:
                    ax.plot(nt.t*1e3,nt.y[i,:],linewidth=linewidth,linestyle=linestyle)
    
        # Plot total number of nuclei as function of time
        if nuclei is True:
            ax.plot(nt.t*1e3,self.totpart(nt.y),'k',linewidth=linewidth,label='Total'+labelapp,linestyle=linestyle)
        else:
            ax.plot(nt.t*1e3,sum(nt.y,axis=0),'k',linewidth=linewidth,label='Total'+labelapp,linestyle=linestyle)
        # Plot initial number of nuclei
        ax.axhline(self.totpart(n0),color='k',linewidth=0.5)
        
        # Unless title is requested, show automated title
        if title is None: # Automated title
            ex=int(log10(ne))
            mu=ne/(10**ex)
            ax.set_title(pretitle+r' $\rm{{T_e}}$={} eV, $\rm{{n_e}}$={}$\times10^{{{}}}$ $\rm{{cm^{{-3}}}}$'.format(Te,mu,ex))
        else: # Requested title
            ax.set_title(title)

        ax.legend(ncol=ncol,loc='best') # Show legend
        ax.set_xlim(0,t*1e3) # Set X-lim in ms
        ax.set_xlabel('Time [ms]') # Show x-label
        ax.set_ylabel(r'Density [$\rm{cm^{-3}}$]') # Show y-label
        # If no ylim is requested set the ylim to the highest plotted line
        if ylim is None:
            # Total particles in this plot
            totmax=max(self.totpart(nt.y))
            # Compare to previous max
            if totmax>ax.get_ylim()[1]:
                # If higher, update
                ax.set_ylim((0,totmax))
            else: # Set bottom to zero
                ax.set_ylim(bottom=0)
        else: # Use requested ylim
            ax.set_ylim(ylim)
        
        # Save to casepath/output/figs/figname.type if figname set
        if savename is not None:
            # Make sure that the direcory exists
            try:
                mkdir('{}/output/figs'.format(self.path))
            except:
                pass
            # Save
            fig.savefig('{}/output/figs/{}.{}'.format(self.path,savename,figtype),dpi=300,edgecolor=None,format=figtype,bbox_inches='tight')


        fig.show() # Show fig
        return fig # Return figure object
            
            




    def plot_uerate(self,ratefile,Te,ne,fig=None,ax=0,figsize=(10,10/1.618033),linestyle='-',ylim=(1e-1,1e7),linewidth=3,labelapp='',savename=None,title=None,pretitle='',figtype='png',color='k',plot='loglog',idx=0,xlim=None,ncol=3):
        ''' Routine plotting the UEDGE rates read from file ratefile
            plot_uerate(ratefile,Te,ne,*keys)
        
            ratefile    -   File with rate data in UEDGE format
            Te          -   Background electron plasma temperature [eV]
            ne          -   Background electron plasma density [cm**-3]
            
            Optional parameters
            plot ('loglgog')-   Plot type from options 'loglog', 'semilogy', 'semilogx', or 'plot'
            idx (0)         -   The index of appearance of the data block in ratefile
            fig (None)      -   Figure object on which to plot. Creates new figure if None or unviable
            ax (0)          -   Which axis of fig to plot on, if fig is viable for plotting
            figsize (10,16) -   Figure size, if created by script
            linestyle ('-') -   Plot line style for the requested spaces
            linewidth (2)   -   Line width of P-space and all lines if N-space; Q-space linewidth=0.5*linewidth
            color ('k')     -   Line color
            title (None)    -   User-defined title, overrides auto-generated title
            pretitle ('')   -   Prefix for auto-generated title (temperature and density of BG plasma)
            labelapp ('')   -   String to append to labels in legend
            ncol (3)        -   Number of columns in the legend
            ylim (1e-1,1e7) -   Requested y-limits
            ylim (None)     -   Requested y-limits,defaults to standard if None
            savename (None) -   Figure savename. Default location path/output/figs. Not saved if None
            figtype ('png') -   Figure type if saved
            
            
        '''
        from matplotlib.pyplot import figure
        from numpy import log10,array
        from os import mkdir
        from CRUM.reactions import REACTION


        rates={}
        # Read the data from ratefile into dict
        datalist=['H0_depl','H0_create','H2_depl','H2_create','H0_ext','H2_ext']
        self.ratedata.read_UE(ratefile,rates,datalist=datalist)
        
        # Create custom reaction
        reactions={}
        for r in rates.keys():
            reactions[r]=REACTION(r,'',[''],[''],rates[r],'UE')
        
        typ,lab='',''
        if idx in [0,1,4]:
            typ+=r'Atomic hydrogen '
            lab+='Atom. H '
        else:
            typ+=r'Molecular hydrogen '
            lab+='Mol. H '
        if idx in [0,2]:
            typ+='depletion'
            lab+='sink'
        elif idx in [1,3]:
            typ+='creation'
            lab+='source'
        else:
            typ+='external source (EIR)'
            lab+='ext. source'


        # Setup plot type
        if isinstance(Te,list):
            x=Te
            y=array([reactions[datalist[idx]].rate(T,None,ne=ne) for T in Te])/(ne**(idx>3))
            xlabel=(r'Electron temperature [eV]')
            ex=int(log10(ne))
            mu=ne/(10**ex) 
            autotitle=pretitle+r' {}, $\rm{{n_e}}$={}$\times10^{{{}}}$ $\rm{{cm^{{-3}}}}$'.format(typ,mu,ex)
        elif isinstance(ne,list): 
            x=ne
            y=array([reactions[datalist[idx]].rate(Te,None,ne=n) for n in ne])
            if idx>3:
                y=y/ne
            xlabel=(r'Electron density [$\rm{cm^{-3}}$]')
            autotitle=pretitle+r' {}, $\rm{{T_e}}$={} [eV]'.format(typ,Te)
        else:
            print('Reaction rate for process "{}" at Te={} and ne={}'.format(datalist[idx],Te,ne))
            return reactions[datalist[idx]].rate(Te,None,ne)
        
        try:
            ax=fig.get_axes()[0]
        except:
            fig=figure(figsize=figsize)
            ax=fig.add_subplot(111)

        if plot=='loglog':
            pl=ax.loglog
        elif plot=='semilogy':
            pl=ax.semilogy
        elif plot=='semilogx':
            pl=ax.semilogx
        elif plot=='plot':
            pl=ax.plot
        else:
            print('Unrecognized plot "{}" requested! Aborting.'.format(plot))


        pl(x,abs(y),linewidth=linewidth,color=color,label=lab+labelapp,linestyle=linestyle)

        if title is None:
            ax.set_title(autotitle)
        else:
            ax.set_title(title)

        ax.legend(ncol=ncol,loc='best')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$\rm{<\sigma v>}$ [$\rm{n_e\cdot V\quad s^{-1}}$]')
        if xlim is not None:
            ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        
        if savename is not None:
            try:
                mkdir('output/figs')
            except:
                pass
            fig.savefig('output/figs/{}.{}'.format(savename,figtype),dpi=300,edgecolor=None,format=figtype,bbox_inches='tight')


        fig.show()
        return fig
