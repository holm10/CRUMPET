# CRUM input parser class readinput.py
# Separated from CRUM.py by holm10
# Changelog
# 200205 - Separated from CRUM.py #holm10
from CRUM.crm import CRM
from CRUM.ratedata import RATE_DATA
from CRUM.tools import TOOLS


# TODO: Object for species?


class CRUMPET(CRM,RATE_DATA,TOOLS):
    '''CRUMPET container object for collsional-radiative models

    This class is used to read and parse the input file, the parts
    which are passed onto the CRM class that sets up the CRM itself.
    
    ...

    Attributes:
    -----------
    ver : srt
        a string defining the code version

    Methods
    -------
    create_UE_rates
    plot_depletionT
    lifetimes
    plot_Et
    plot_nt
    get_UE_reactions
    plot_UE_nrate
    calc_UE_Erate
    plot_UE_Erate
    spectrum
    plotrate
    steady_state
    

    
    species : dict
        Each entry is a dictionary containing information on the 
        initial density and potential level of the CRM  species. The 
        dictionary keys are used to identify the species
    bg : dict
        As for 'species', but for the plasma background
    NP : int (default 2)
        Number of P-species (time-dependently evaluated species) to be
        considered in the sumulation. First NP species defined in the
        input deck is considered.
    vmax : int
        Number of maxium vibrational states inculded in the CRM
    nmax : int
        Number of maximum atomic electronic states included in the CRM
    verbose : bool
        Switch whether to show verbose output in the terminal. The 
        same information os written to the log
    path :  str (default '.')
        Path to the parent folder relative to which subsequent paths
        are defined
    '''


    def __init__(self,  fname='input/CRUM.dat',path='.',vmax=14,nmax=8,
                        verbose=False,NP=2):
        '''
        Parameters
        ----------
        fname : str, optional (default: 'input/CRUM.dat')
            path to CRUM input file relative to parameter path
        path : str, optional (default: '.')
            path to which remaing paths are relative
        vmax : int, optional (default: 14)
            number of excited vibrational molecular states considered
        nmax : int, optional (default: 8)
            number of electronic excited states considered
        verbose : boolean, optional (default: False)
            displays progress to prompt (written to a log-file too)
        NP : int, optional (default: 2)
            P-space size, chosen as the first NP entries of the species
            card
        '''
        from numpy import zeros,pi
        from os import getcwd

        self.ver='V0.2'



        # Read the input file into a list
        lines,cards,subcards=self.file2list(path,fname)

                
        # %%%%%%%%%% LOOP THROUGH THE DEFINED CARDS %%%%%%%%%%%
        for i in range(len(cards)-1):


            # %%%% Store species from species card to temporary list %%%%
            if lines[cards[i]].split()[1].upper()=='SPECIES':
                species={}

                # Loop through each subcard
                for j in range( subcards.index(cards[i])+1,
                                subcards.index(cards[i+1])  ):
                    buff={}
                    ind=subcards[j]

                    
                    # Loop through the lines following the  declaration 
                    # looking for the potential energy
                    for l in range(subcards[j]+1,subcards[j+1]): 
                        buff[lines[l].split()[0].strip()] = float(
                                                lines[l].split()[1].strip() )


                    species[lines[ind].split()[1].strip()]=buff
                    
            # %%% Store the reactions requested to temporary reaction list %%%
            elif lines[cards[i]].split()[1].upper()=='REACTIONS':
                reactionlist=[]

                # Look between this card and the next: two lines per reactions
                for k in range( subcards.index(cards[i])+1,
                                subcards.index(cards[i+1])  ):

                    # Store the reaction name, discard subcard marker
                    rname=lines[subcards[k]].split()[1]   

                    # Check if custom reaction deck called
                    if rname.upper()=='CUSTOM': 
                        # Store name and path to custom file
                        reactionlist.append( [  rname, 
                                                lines[subcards[k]+1].strip()
                                             ]   ) 

                    
                    # Regular reaction
                    else:
                        # Define energy sink/source
                        if subcards[k]+2==subcards[k+1]: K='0' # No sink/source

                        for l in range( subcards[k]+2,
                                        subcards[k+1]   ): # Read sink/source

                            K=lines[l].strip().split('=')[-1] 
                        
                        # Extract the reactants and fragments from the input
                        reactant,fragment = [ x.strip().split(' + ') for x 
                                        in lines[subcards[k]+1].split(' > ') ]


                        reactionlist.append( [  rname, reactant, fragment, K ])
                        


            # %%%%%% Set paths of the standard rate data files %%%%%
            elif lines[cards[i]].split()[1].upper()=='RATES':
                for j in range(cards[i]+1,cards[i+1]): # Look between this card and the next
                    # Extract the database type and path
                    [db,nm]=lines[j].split()
                    db=db.upper()
                    if db=='ADAS':      ADAS=   nm
                    elif db=='UE':      UE=     nm
                    elif db=='AMJUEL':  amjuel= nm
                    elif db=='HYDHEL':  hydhel= nm
                    elif db=='H2VIBR':  h2vibr= nm
                    else:
                        print('''   Unrecognized database type {}. 
                                    Ignoring.'''.format(db))
                        continue

                # Create RATE_DATA class object based on the above files paths
                RATE_DATA.__init__(self,amjuel,hydhel,h2vibr,ADAS,UE,path)



            # %%%% Define CRM settings %%%%
            elif lines[cards[i]].split()[1].upper()=='SETTINGS':

                # Search this card only
                for j in range(cards[i]+1,cards[i+1]): 
                    # Create setting based on card type
                    try:
                        [c,setting,value]=lines[j].upper().split()
                    except:
                        [c,setting]=lines[j].upper().split()

                    if c!='*': continue # Ensure that we are reading subcards
                    elif setting=='VMAX':       vmax=int(value)
                    elif setting=='NMAX':       nmax=int(value)
                    elif setting=='VERBOSE':    verbose=bool(int(value))
                    elif setting=='NP':         Np=int(value)
                    elif setting=='N0':
                        # Loop through all defined initial densities
                        for k in range( j+1, subcards[subcards.index(j)+1 ]): 

                            # Store the density in the location of n0 that 
                            # corresponds to the species requested
                            species[lines[k].split()[0]]['n']=float(
                                                        lines[k].split()[1]) 

                    else:
                        print(   '''Unrecognized setting {}. 
                                    Ignoring.'''.format(setting))
                        continue

            # %%%% Define the plasma background species %%%%
            elif lines[cards[i]].split()[1].upper()=='BACKGROUND':
                bg={}
                # Loop through each subcard
                for j in range( subcards.index(cards[i])+1,
                                subcards.index(cards[i+1])  ):
                    buff={}
                    ind=subcards[j]

                    # Loop through the lines following the declaration 
                    # looking for the background species potential
                    for l in range( subcards[j]+1, subcards[j+1]): 
                        buff[lines[l].split()[0].strip().upper()]=float(
                                            lines[l].split()[1].strip())


                    bg[lines[ind].split()[1].strip()]=buff

            # %%%%% Unknown card, abort %%%%%%
            else:
                print( '''Unknown card "{}": 
                            aborting!'''.format(lines[cards[i]].split()[1]))
                return



        # %%%% Set up the CRM %%%%
        CRM.__init__(self,  species,bg,reactionlist,verbose,Np,
                            path,vmax,nmax,self.reactions           )
     






    def create_UE_rates(self,fname='ue',E=0.1,Sext=True,h0h2=['H(n=1)','H2(v=0)'],Tm=False,Ton=False,rad=False):
        ''' Script that writes UEDGE rates to self.path/fname.dat
            create_UE_rates(fname='uerates')

            Optional parameters
            E (0.1)     -   target particle energy [eV]
            Sext (True) -   Include external source (from background plasma reactions into CRM species)
            h0h2 (['H(n=1)','H2(v=0)']) - List of the H0 and H2 at their electronic and vibrational ground state handles used in the input
            '''
        from numpy import zeros,sum,array
        from os import mkdir
        from datetime import datetime
        from getpass import getuser

        # Ensure that the CRM format coincides with the assumed format
        if list(self.species)[:2]!=h0h2:
            print('This model assumes the first species to be H0 and the second to be H2 in their electronic and vibrational ground states.\nThe first input species do not match {} and {}, respectively. Aborting.'.format(h0h2[0],h0h2[1]))
            return

        if self.Np!=2:
            print('Presently, UEDGE only considers two neutral hydrogenic species, H0 and H2. The P-space of the model is larger than thisand, thus, is incompatible with the UEDGE model. No rates can be written for Np!=2. Aborting.')
            return

        # Create deinsty and temperature points in log-log space as in existing UEDGGE rate files
        ret=zeros((15,60,self.Np,self.Np))
        retE=zeros((15,60,5,self.Np))
        ext=zeros((15,60,self.Np))
        extE=zeros((15,60,5))
        # Calculate the Greenland (Np) space rates for the T-n space
        print('Creating UEDGE rate output')
        for i in range(15): # 15 density points in log-log space
            print('    Density step {}/15'.format(i+1))
            for j in range(60): # 60 temperature points in log-log space
                print('        Temperature step {}/60'.format(j+1))
                Te,ne=10**(-1.2+j/10),10**(10+0.5*i)
                ret[i,j,:,:],ext[i,j,:],_=self.gl_crm(*self.M(Te,ne,Te,ne,E,write=False),Sext=Sext) # Store to matrix
                U=self.Sgl(Te,ne,Te,ne,E,rad,Tm,write=False,Ton=Ton) # Don't include erl1/erl2 radiation in the rates as these are handled by UEDGE
                # Include T-losses or not?
                retE[i,j,:,:]=array([ sum(U[0][0],axis=0), sum(U[1][0],axis=0), sum(U[2][0],axis=0), sum(U[3][0],axis=0), sum(U[4][0],axis=0)])
                #print(len(retE.shape))
                extE[i,j,:]=array([ sum(U[0][1],axis=0), sum(U[1][1],axis=0), sum(U[2][1],axis=0), sum(U[3][1],axis=0), sum(U[4][1],axis=0) ])
        

        # Ensure output directory exists 
        try:
            mkdir('{}/output'.format(self.path))
        except:
            pass
        
        # Write UEDGE data
        print(' Writing UEDGE reaction rates to {}'.format(fname+'_nrates'))
        with open('{}/output/{}.dat'.format(self.path,fname+'_nrates'),'w') as f:
            # Headers and positions in Greenland matrices
            setups= [   [' H0 depl. Rate(jt,jn) (s**-1)  Te(jt) = 10**(-1.2 + (jt-1)/10) jt= 1,60\n', ret[:,:,0,0]    ],
                        [' H2->H0 Rate(jt,jn) (s**-1)\n', ret[:,:,0,1]                                                ],
                        [' H2 depl. Rate(jt,jn) (s**-1)\n', ret[:,:,1,1]                                              ],
                        [' H2 creation Rate(jt,jn) (s**-1)\n', ret[:,:,1,0]                                           ],
                        [' H0 external source Rate(jt,jn) (s**-1)\n', ext[:,:,0]                                     ],
                        [' H2 external source Rate(jt,jn) (s**-1)\n', ext[:,:,1]                                     ]
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
            f.write('\n'+datetime.now().strftime("%Y-%m-%d %H:%M:%S")+'by user '+getuser()+' using CRUMPET '+self.ver)


        print(' Writing UEDGE energy rates to {}'.format(fname+'_Erates'))
        with open('{}/output/{}.dat'.format(self.path,fname+'_Erates'),'w') as f:
            # Headers and positions in Greenland matrices
            setups= [   [' e-loss (H) (eV s**-1)  Te(jt) = 10**(-1.2 + (jt-1)/10) jt= 1,60\n', retE[:,:,0,0]    ],
                        [' e-loss (H2) (eV s**-1)\n', retE[:,:,0,1]                                                ],
                        [' ext. e-loss (eV s**-1)\n', extE[:,:,0]                                                ],
                        [' ia source (H) (eV s**-1)\n', retE[:,:,1,0]                                                ],
                        [' ia source (H2) (eV s**-1)\n', retE[:,:,1,1]                                                ],
                        [' ext ia source (eV s**-1)\n', extE[:,:,1]                                                ],
                        [' Epot source (H) (eV s**-1)\n', retE[:,:,2,0]                                                ],
                        [' Epot source (H2) (eV s**-1)\n', retE[:,:,2,1]                                                ],
                        [' ext Epot source (eV s**-1)\n', extE[:,:,2]                                                ],
                        [' rad source,a (H) (eV s**-1)\n', retE[:,:,3,0]                                                ],
                        [' rad source,a (H2) (eV s**-1)\n', retE[:,:,3,1]                                                ],
                        [' ext rad source,a (eV s**-1)\n', extE[:,:,3]                                                ],
                        [' rad source,m (H) (eV s**-1)\n', retE[:,:,4,0]                                                ],
                        [' rad source,m (H2) (eV s**-1)\n', retE[:,:,4,1]                                                ],
                        [' ext rad source,m (eV s**-1)\n', extE[:,:,4]                                                ],
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
            f.write('\n'+datetime.now().strftime("%Y-%m-%d %H:%M:%S")+'by user '+getuser()+' using CRUMPET '+self.ver)


    def plot_depletionT(self,key,n,Te,Ti=None,color='r', fig=None,figsize=(10,6.25),E=0.1,style=['-','-.','--',':']):
        ''' T-dependent plot of the total rate for species '''
        from matplotlib.pyplot import figure
        if Ti is None: Ti=Te

        try:
            ax=fig.get_axes()[0]
        except:
            fig=figure(figsize=figsize)
            ax=fig.add_subplot(111)

        try:
            len(n)
            nl=n
        except:
            nl=[]
            nl.append(n)
            

        for n in range(len(nl)):
            y=[]
            for i in range(len(Te)):
                y.append(self.species_rate(key,Te[i],nl[n],Ti[i],n,E))
            ax.loglog(Te,y,color=color,linestyle=style[n])
                


        return fig 


    def lifetimes(self,Te,ne,E=0.1,Ti=None,ni=None,Sext=True,n=None,reassociation=1e3):
        ''' Returns and prints a list of radiative lifetimes for each species '''
        from numpy import log,ones,matmul,polyfit,zeros,exp,sum
        from scipy.optimize import fsolve,minimize
        from scipy.linalg import norm
        from scipy.integrate import solve_ivp
        from matplotlib.pyplot import figure


        def dt(t,n,mat,ext,n0):
            #ntot=self.totpart(n)
            #dn=self.totpart(n)-self.totpart(n0)
            #ext[1]=-dn*fac
            #print(n[1],dn)
            #fac=0*1e3#-1e-2
            '''ext[0]=-(matmul(mat,n))[0]'''
            
            #ext[1]=-0.5*ext[0]
            #ext[1]=-ext[0]*0.1
            #ext[0]=-fac*n[0]
            #ext[1]=-0.99*self.totpart(ext)
            #print(self.totpart(ext))
            #ext[1]=-0.5*ext[0]
            #print(self.totpart(n),sum(ext))
            return matmul(mat,n)+ext

        def ss(t,Y):
            global ss_n 
            if norm(Y-ss_n)<7e5:
                return 0
            else:
                return 1
        ss.terminal=True

        def findfac(fac,n,mat,ext,n0):
            nend_sum=self.totpart(solve_ivp(lambda x,y: dt(x,y,mat,ext,n,fac),(0,1),n,'LSODA',dense_output=True).y[:,-1])
            
            return abs(self.totpart(n0)-nend_sum)
      



        na=2.22793E+16
        nm=7.64786E+17
        diva=-8.194620066655537e+18
        divm=-1.034434079121622e+18
        srca=0.00000E+00
        srcm=6.24142E+18
        bga=2.42937E+08
        bgm=9.19138E+08
        psor=-4.00924E+16
        psorgc=2.63178E+12
        vol=1e2

          
        na=5.03524E+15
        nm=1.20777E+17
        diva=-1.01910E+15
        divm=-2.10175E+11
        srca=0.00000E+00
        srcm=6.24142E+18
        bga=2.31957E+08
        bgm=5.15052E+09
        psor=-8.65157E+15
        psorgc=1.92209E+15
        vol=1e4

        na=3.42356E+15
        nm=3.62223E+17
        diva=-3.69215E+18
        divm=-4.89929E+18
        srca=0.00000E+00
        srcm=6.24142E+18
        bga=2.15803E+04
        bgm=5.00209E+08
        psor=-5.47270E+11
        psorgc=6.16478E+15
        vol=1e4

        na=1.48518E+16
        nm=3.59922E+18
        diva=-2.81098E+18
        divm=-4.86823E+18
        srca=0.00000E+00
        srcm=6.24142E+18
        bga=2.31958E+06
        bgm=5.15054E+07
        psor=-2.55184E+14
        psorgc=1.92210E+13
        vol=1e2

        na=3.23058E+15
        nm=1.81807E+16
        diva=-1.10649E+19
        divm=-2.45590E+11
        srca=0.00000E+00
        srcm=6.24142E+18
        bga=8.34627E+09
        bgm=4.63451E+10
        psor=-1.99728E+17
        psorgc=6.25988E+14
        vol=1e4

        # TODO: Check that volumes go right: CM to CM
        # TODO: Implement approximation of pumping 

        if n is None: n=self.n0() # Use input n0 as default unless explict n0 requested
        #NN=len(self.species)*len(self.species)
        n[0]=1e-6*na
        n[1]=1e-6*nm

       
        mat,ext=self.M(Te,ne,Ti,ni,E,write=False) # Get the full rate matrix

        # TODO: use Greenland CRM to assess the SS - the same assumptions
        ''' Assume re-association to achieve SS '''
        #mat[0,0]=-reassociation
        #mat[1,0]=-mat[0,0]*0.5 # Conserve nucleii
        ext[0]=(psorgc+diva+srca+bga)
        ext[1]=(divm+srcm+bgm)
        mat[0,0]=(psor/ne)


        M,T,D=self.gl_crm(mat,ext,matrices=True) # Get matrices and eigenvalues/-vectors

        Meff,extp,n0p=self.gl_crm(mat,ext,n=n) # Get matrices and eigenvalues/-vectors
        #print(matmul(mat[:,1:],n[1:])[0],mat[1,1],Meff[1,1],n)
        for i in D:
            print('{:.3E}'.format(-1/i))
        #fac=minimize(findfac,5e4,(n,mat,ext,n))
        
        # Supply molecules to achieve SS
        '''ext[1]=sum(n)/0.2'''
        # Simulate to SS
        #ss_sol=solve_ivp(lambda x,y: dt(x,y,Meff,extp,n0p),(0,0.1),n0p,'LSODA',dense_output=True)
        ss_sol=solve_ivp(lambda x,y: dt(x,y,M,ext,n),(0,10),n,'LSODA',dense_output=True)
        f=figure(figsize=(7,7))
        ax=f.add_subplot(111)
        for i in range(len(n)):
            line='--'
            if 'H2' in self.slist[i]: line='-'

            ax.semilogx(ss_sol.t,ss_sol.y[i,:],line,label=self.slist[i]) 
            if i==0:
                ax.semilogx(ss_sol.t,ss_sol.y[i,:],'k.',label=self.slist[i]) 
        ax.legend(ncol=3)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Density [cm**-3]')
        #print(matmul(mat[:,1:],n[1:])[0],mat[0,1],Meff[0,1],mat[1,1],Meff[1,1],n[1])
        print('UEDGE molecular content={:.5E} m**-3'.format(nm))
        print('CRUMPET molecular content={:.5E} m**-3'.format(1e6*self.totmol(ss_sol.y[:,-1])))
        return
        
        global ss_n
        ss_n=ss_sol.y[:,-1]

        ret=[1e20]
        # Perturb each species
        for i in range(1,len(n)):
            try:
                pert=zeros((len(n),))
                pert[i]=1e6
                sol=solve_ivp(lambda x,y: dt(x,y,mat,ext,n),(0,1e-3),ss_n+pert,'LSODA',dense_output=True,events=ss)
                a,b=polyfit(sol.t,log(sol.y[i,]),1)
                print('tau({})'.format(self.slist[i]).ljust(20,' ')+'= {:.3E}'.format(-1/a))
            except: 
                pert=zeros((len(n),))
                pert[i]=1e6
                sol=solve_ivp(lambda x,y: dt(x,y,mat,ext,n),(0,1e-3),ss_n+pert,'LSODA',dense_output=True,events=ss)
                a,b=polyfit(sol.t,log(sol.y[i,]),1)
                print('tau({})'.format(self.slist[i]).ljust(20,' ')+'= {:.3E}'.format(-1/a))

            ret.append(-1/a)
                #sol=solve_ivp(lambda x,y: dt(x,y,mat,ext),(0,1e-3),ss_n+pert,'LSODA',dense_output=True)
                #for i in range(len(n)):
                #    ax.plot(sol.t,sol.y[i,:],label=self.slist[i]) 
        #f=figure(figsize=(7,7))
        #ax=f.add_subplot(111)
        #for i in range(len(n)):
        #    ax.plot(sol.t,sol.y[i,:],label=self.slist[i]) 
        #ax.plot(sol.t,sol.y[1,0]*exp(sol.t*a),'k-')
        #ax.plot(sol.t,sol.y[1,:],'.-')     
        #ax.plot(sol.t,(sol.t*a+b),'k-')
        #ax.plot(sol.t,log(sol.y[1,:]),'.-')     
        #ax.legend(ncol=3)
        print('====')

        for i in range(1,len(n)):
            print('tau({})'.format(self.slist[i]).ljust(20,' ')+'= {:.3E}'.format(-1/mat[i,i]))

        return ret 


    def plot_Et(self,t,Te,ne,E=0.1,Ti=None,ni=None,Sext=True,Np=True,Nq=False,N=False,fig=None,n0=None,ax=0,figsize=(7+7,7*1.618033),linestyle='-',ylim=None,linewidth=2,labelapp='',qlabel=False,savename=None,title=None,pretitle='',figtype='png',color=None,nuclei=True,Tm=False,rad=True,gl=False,n=None,suptitle=None,ylimE=None,Ton=True):
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
        from numpy import log10,sum,matmul,zeros,linspace,amax,where
        from numpy.linalg import inv
        from os import mkdir
       
#        for r in self.reactions: 
#            print('Reaction {}: S_r={}, S_g={}, S_V={}, S_e={}'.format(r.name,r.S_r,r.S_g,r.S_V,r.S_e))
        # Get axis handle or create figure
        try:
            ax=fig.get_axes()[0]
        except:
            fig=figure(figsize=figsize)
        
        ev=1.602e-19
        # If no initial distribution is requested, use the one specified in the input
        if n0 is None:
            n0=self.n0()

        Neq=False
        if gl is True:
            Neq=Nq
            Np,Nq,N=True,False,False

        CRM=self.dEdt(t,Te,ne,Ti,ni,E,Tm,rad,Sext,gl=gl,n=n,Qres=((Nq is True) or (N is True) or (Neq is True)),Ton=Ton) # Set up Greenland model 

        if color is None: # Define color sequence up to 10 unless specific sequence requested
            color=[ 'b', 'r', 'm', 'c', 'darkgreen', 'gold', 'brown' ,'lime', 'grey', 'orange' ]


        t=linspace(0,t,800)


        if (Nq is False) and (N is False) and (Neq is False):
            ''' Only plot P-species and contributions, don't resolve energy '''


            ax=fig.add_subplot(311)

            part=CRM.sol(t)[5:]

            # Plot P-space species if requested
            if Np:
               for i in range(self.Np):
                    ax.plot(t*1e3,part[i,:],linewidth=linewidth,color=color[i],label=list(self.species)[i]+labelapp,linestyle=linestyle)
           # Plot total number of nuclei as function of time
            if nuclei is True:
                ax.plot(t*1e3,self.totparticles(part),'k',linewidth=linewidth,label='Total'+labelapp,linestyle=linestyle)
            else:
                ax.plot(t*1e3,sum(part,axis=0),'k',linewidth=linewidth,label='Total'+labelapp,linestyle=linestyle)
            # Plot initial number of nuclei
            ax.axhline(self.totparticles(n0),color='k',linewidth=0.5)
            

            ax.legend(ncol=self.Np+1,loc='best') # Show legend
            ax.set_xlim(0,t[-1]*1e3) # Set X-lim in ms
            ax.set_ylabel(r'Density [$\rm{cm^{-3}}$]') # Show y-label
            # If no ylim is requested set the ylim to the highest plotted line
            if ylim is None:
                # Total particles in this plot
                totmax=max(self.totparticles(CRM.sol(t)[5+len(self.n0())*(gl==True):]))
                # Compare to previous max
                if totmax>ax.get_ylim()[1]:
                    # If higher, update
                    ax.set_ylim((0,totmax))
                else: # Set bottom to zero
                    ax.set_ylim(bottom=0)
            else: # Use requested ylim
                ax.set_ylim(ylim)


            ax=fig.add_subplot(312)

            E=[-CRM.sol(t)[0],CRM.sol(t)[1],CRM.sol(t)[3],CRM.sol(t)[4],CRM.sol(t)[2]]
            Elabel=[r'$\mathrm{S_{e-loss}}$',r'$\mathrm{S_{ia}}$',r'$\mathrm{S_{rad,a}}$',r'$\mathrm{S_{rad,m}}$',r'$\mathrm{S_{pot}}$']

            for j in range(5):
                ax.plot(t*1e3,E[j]*ev,linewidth=linewidth,color=color[j],label=Elabel[j]+labelapp,linestyle=linestyle)

            ax.set_ylabel(r'Power [W]') # Show y-label
            ax.legend(ncol=4,loc='best') # Show legend
            ax.axhline(y=0,color='k')
            ax.set_xlim(0,t[-1]*1e3) # Set X-lim in ms       '''
            if ylimE is not None:
                ax.set_ylim(ylimE)

            # Show balances
            ax=fig.add_subplot(313)
            exp=round(log10(max(abs(sum(CRM.sol(t)[0:5],axis=0))*ev)))
            ax.plot(t*1e3,sum(CRM.sol(t)[0:5],axis=0)*ev*(10**(-exp)),linewidth=linewidth,color='r',label=r'$\mathrm{\Delta E}$'+r' [$10^{{{:d}}}$W]'.format(int(exp))+labelapp,linestyle=linestyle)
            ax.set_xlabel('Time [ms]') # Show x-label
            ax.set_ylabel(r'Conservation') # Show y-label
            ax.legend(loc='best') # Show legend
            ax.set_xlim(0,t[-1]*1e3) # Set X-lim in ms

        else:
            Nt=len(self.species)

            # Particle balance        
            ax=fig.add_subplot(421)

            part=CRM.sol(t)[-Nt:]
            if Neq is True:
                part=CRM.sol(t)[-self.Np:]

            # Plot P-space species 
            for i in range(self.Np):
                ax.plot(t*1e3,part[i,:],linewidth=linewidth,color=color[i],label=list(self.species)[i]+labelapp,linestyle=linestyle)
            # Plot Q-space species if requested
            if Nq:
                for i in range(self.Np,len(self.species)):
                    if qlabel:
                        ax.plot(t*1e3,part[i],linewidth=linewidth*(0.5)**Np,label=list(self.species)[i]+labelapp,linestyle=linestyle)
                    else:
                        ax.plot(t*1e3,part[i],linewidth=linewidth*(0.5)**Np,linestyle=linestyle)
            # Plot all species w/ automated colors if requested
            if N:
               for i in range(len(self.species)):
                    if qlabel:
                        ax.plot(t*1e3,CRM.sol(t)[i,:],linewidth=linewidth,label=list(self.species)[i-4]+labelapp,linestyle=linestyle)
                    else:
                        ax.plot(t*1e3,CRM.sol(t)[i,:],linewidth=linewidth,linestyle=linestyle)
        
            # Plot total number of nuclei as function of time
            if nuclei is True:
                ax.plot(t*1e3,self.totparticles(part),'k',linewidth=linewidth,label='Total'+labelapp,linestyle=linestyle)
            else:
                ax.plot(t*1e3,sum(part,axis=0),'k',linewidth=linewidth,label='Total'+labelapp,linestyle=linestyle)
            # Plot initial number of nuclei
            ax.axhline(self.totparticles(n0),color='k',linewidth=0.5)

            ax.legend(loc='best') # Show legend
            ax.set_xlim(0,t[-1]*1e3) # Set X-lim in ms
            ax.set_xlabel('Time [ms]') # Show x-label
            ax.set_ylabel(r'Density [$\rm{cm^{-3}}$]') # Show y-label
            ax.set_title('Species density')
            # If no ylim is requested set the ylim to the highest plotted line
            if ylim is None:
                # Total particles in this plot
                totmax=max(self.totparticles(part))
                # Compare to previous max
                if totmax>ax.get_ylim()[1]:
                    # If higher, update
                    ax.set_ylim((0,totmax))
                else: # Set bottom to zero
                    ax.set_ylim(bottom=0)
            else: # Use requested ylim
                ax.set_ylim(ylim)



            
            # Energy 
            Etitle=['Electron loss','Ion/atom source','Potential','Radiation,a','Radiation,m']

            for e in range(5):
                E=((-1)**(e==0))*CRM.sol(t)[Nt*e:Nt*(e+1)]*ev
                
                
                ax=fig.add_subplot(4,2,2+e) 

                # Plot P-space species if requested
                if Np:
                   for i in range(self.Np):
                        if abs(E[i].max())!=0:
                            ax.plot(t*1e3,E[i],linewidth=linewidth,color=color[i],label=list(self.species)[i]+labelapp,linestyle=linestyle)
                # Plot Q-space species if requested
                if Nq or Neq:
                    for i in range(self.Np,len(self.species)):
                        if abs(E[i].max())!=0:
                            ax.plot(t*1e3,E[i],linewidth=linewidth*(0.5)**Np,label=list(self.species)[i]+labelapp,linestyle=linestyle)
                # Plot all species w/ automated colors if requested
                if N:
                   for i in range(len(self.species)):
                        if abs(E[i].max())!=0:
                            ax.plot(t*1e3,E[i],linewidth=linewidth,label=list(self.species)[i]+labelapp,linestyle=linestyle)
                Etot=sum(E,axis=0)
                ax.plot(t*1e3,Etot,'k',linewidth=linewidth,label='Total'+labelapp,linestyle=linestyle)
                ax.set_title(Etitle[e])
                ax.set_ylabel(r'Power [W]') # Show y-label
                ax.set_xlim(0,t[-1]*1e3) # Set X-lim in ms       '''
                ax.set_xlabel('Time [ms]') # Show x-label
                if e==0: 
                    Emax=max(Etot)
                    axmax=ax.get_ylim()[1]
                if ylimE is not None:
                    ax.set_ylim(ylimE)
                else:
                    if Emax>axmax:
                        # If higher, update
                        ax.set_ylim(top=Emax)

                box = ax.get_position()
                ax.set_position([box.x0, box.y0 - box.height * 0,box.width, box.height * 0.85])
                if e==4:
                    ax.legend(loc='lower center', bbox_to_anchor=(-0.1, -0.9),frameon=False, ncol=9)

            # Show balances
            ax=fig.add_subplot(414)
        
            exp=round(log10(max(abs(sum(CRM.sol(t)[:5*Nt],axis=0))*ev)))
            ax.plot(t*1e3,sum(CRM.sol(t)[0:5*Nt],axis=0)*ev*(10**(-exp)),linewidth=linewidth,color='r',label=r'$\mathrm{\Delta E}$'+r' [$10^{{{:d}}}$W]'.format(int(exp))+labelapp,linestyle=linestyle)
            ax.set_xlabel('Time [ms]') # Show x-label
            ax.set_ylabel(r'Conservation') # Show y-label
            ax.legend(loc='best') # Show legend
            ax.set_xlim(0,t[-1]*1e3) # Set X-lim in ms
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 - box.height * 0.3,box.width, box.height * 0.7])




        # Unless title is requested, show automated title
        if suptitle is None: # Automated title
            ex=int(log10(ne))
            mu=ne/(10**ex)
            fig.suptitle(pretitle+r' $\rm{{T_e}}$={} eV, $\rm{{n_e}}$={}$\times10^{{{}}}$ $\rm{{cm^{{-3}}}}$'.format(Te,mu,ex))
        else: # Requested title
            fig.suptitle(title)






        
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
            n0=self.n0()

        # Check what model to use, Greenland or full
        if gl is True: # Greenland
            nt=self.gl_nt(Te,ne,t,Ti,ni,E,n0,Sext) # Solve ODE
            Np,Nq,N=True,False,False # Np only option
        else: # Full model
            nt=self.full_nt(Te,ne,t,Ti,ni,E,n0,Sext) # Solve ODE

        if color is None: # Define color sequence up to 10 unless specific sequence requested
            color=[ 'b', 'r', 'm', 'c', 'darkgreen', 'gold', 'brown' ,'lime', 'grey', 'orange' ]

        # Plot P-space species if requested
        if Np:
           for i in range(self.Np):
                ax.plot(nt.t*1e3,nt.y[i,:],linewidth=linewidth,color=color[i],label=list(self.species)[i]+labelapp,linestyle=linestyle)
        # Plot Q-space species if requested
        if Nq:
            for i in range(self.Np+1,len(self.species)):
                if qlabel:
                    ax.plot(nt.t*1e3,nt.y[i,:],linewidth=linewidth*(0.5)**Np,label=list(self.species)[i]+labelapp,linestyle=linestyle)
                else:
                    ax.plot(nt.t*1e3,nt.y[i,:],linewidth=linewidth*(0.5)**Np,linestyle=linestyle)
        # Plot all species w/ automated colors if requested
        if N:
           for i in range(len(self.species)):
                if qlabel:
                    ax.plot(nt.t*1e3,nt.y[i,:],linewidth=linewidth,label=list(self.species)[i]+labelapp,linestyle=linestyle)
                else:
                    ax.plot(nt.t*1e3,nt.y[i,:],linewidth=linewidth,linestyle=linestyle)
    
        # Plot total number of nuclei as function of time
        if nuclei is True:
            ax.plot(nt.t*1e3,self.totparticles(nt.y),'k',linewidth=linewidth,label='Total'+labelapp,linestyle=linestyle)
        else:
            ax.plot(nt.t*1e3,sum(nt.y,axis=0),'k',linewidth=linewidth,label='Total'+labelapp,linestyle=linestyle)
        # Plot initial number of nuclei
        ax.axhline(self.totparticles(n0),color='k',linewidth=0.5)
        
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
            totmax=max(self.totparticles(nt.y))
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
            
            

    def get_UE_reactions(self,nratefile,eratefile):
        from CRUM.reactions import REACTION
        rates,reactions={},{}
        # Read the data from ratefile into dict
        datalist=['H0_depl','H0_create','H2_depl','H2_create','H0_ext','H2_ext']
        self.read_UE(nratefile,rates,datalist=datalist,path=self.path)
        # Create custom reaction
        reactions={}
        for r in rates.keys():
            reactions[r]=REACTION(r,'',rates[r],'UE',['',0,'',0,[0,0,0,0]])

        # Read the data from ratefile into dict
        datalist=['Hel','H2el','extel','Hia','H2ia','extia','HV','H2V','extV','Hga','H2ga','extga','Hgm','H2gm','extgm',]
        try:
            self.read_UE(eratefile,rates,datalist=datalist,path=self.path)

            # Create custom reaction
            for r in rates.keys():
                reactions[r]=REACTION(r,'',rates[r],'UE',['',0,'',0,[0,0,0,0]])
        except:
            pass

        return reactions






    def plot_UE_nrate(self,ratefile,Te,ne,fig=None,ax=0,figsize=(10,10/1.618033),linestyle='-',ylim=(1e-1,1e7),linewidth=3,labelapp='',savename=None,title=None,pretitle='',figtype='png',color='k',plot='loglog',idx=0,xlim=None,ncol=3,fac=1,E=0.1,Ton=False,rad=True,Tm=False):
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
        from numpy import log10,array,ndarray
        from os import mkdir
        from CRUM.reactions import REACTION


        try:
            rates={}
            # Read the data from ratefile into dict
            datalist=['H0_depl','H0_create','H2_depl','H2_create','H0_ext','H2_ext']
            self.read_UE(ratefile,rates,datalist=datalist)
            # Create custom reaction
            reactions={}
            for r in rates.keys():
                reactions[r]=REACTION(r,'',rates[r],'UE',['',0,'',0,[0,0,0,0]])
        except:
            reactions=self.calc_UE_nrate(Te,ne,E,Ton,rad,Tm)
        
        
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

        if isinstance(Te,list) or isinstance(Te,ndarray):
            x=Te
            try:
                y=array([reactions[datalist[idx]].rate(T,None,ne=ne) for T in Te])/(ne**(idx>3))
            except:
                y=reactions[datalist[idx]]/(ne**idx>3)
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


        pl(x,fac*abs(y),linewidth=linewidth,color=color,label=lab+labelapp,linestyle=linestyle)

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

    def calc_UE_nrate(self,Te,ne,E,Ton,rad,Tm):
        from numpy import zeros,sum,array
        datalist=['H0_depl','H0_create','H2_depl','H2_create','H0_ext','H2_ext']
        ret=zeros((len(Te),self.Np,self.Np))
        ext=zeros((len(Te),self.Np))
        for i in range(len(Te)):
            T=Te[i]
            ret[i,:,:],ext[i,:],_=self.gl_crm(*self.M(T,ne,T,ne,E,write=False),Sext=True) # Store to matrix

        setups= [   ['H0_depl', ret[:,0,0]    ],
                    ['H0_create', ret[:,0,1]                                                ],
                    ['H2_depl', ret[:,1,1]                                              ],
                    ['H2_create', ret[:,1,0]                                           ],
                    ['H0_ext', ext[:,0]                                     ],
                    ['H2_ext', ext[:,1]                                     ]
                ]
        reactions={}
        for l in setups: # Loop through all species contributions
            reactions[l[0]]=l[1]
        return reactions

    def calc_UE_Erate(self,Te,ne,E,Ton,rad,Tm):
        from numpy import zeros,sum,array
        datalist=['Hel','H2el','extel','Hia','H2ia','extia','HV','H2V','extV','Hga','H2ga','extga','Hgm','H2gm','extgm',]
        retE=zeros((len(Te),5,self.Np))
        extE=zeros((len(Te),5))
        for i in range(len(Te)):
            T=Te[i]
            U=self.Sgl(T,ne,T,ne,E,rad,Tm,write=False,Ton=Ton) # Don't include erl1/erl2 radiation in the rates as these are handled by UEDGE
            # Include T-losses or not?
            retE[i,:,:]=array([ sum(U[0][0],axis=0), sum(U[1][0],axis=0), sum(U[2][0],axis=0), sum(U[3][0],axis=0), sum(U[4][0],axis=0)])
            #print(len(retE.shape))
            extE[i,:]=array([ sum(U[0][1],axis=0), sum(U[1][1],axis=0), sum(U[2][1],axis=0), sum(U[3][1],axis=0), sum(U[4][1],axis=0) ])


        setups= [   ['Hel', retE[:,0,0]    ],
                    ['H2el', retE[:,0,1]                                                ],
                    ['extel', extE[:,0]                                                ],
                    ['Hia', retE[:,1,0]                                                ],
                    ['H2ia', retE[:,1,1]                                                ],
                    ['extia', extE[:,1]                                                ],
                    ['HV', retE[:,2,0]                                                ],
                    ['H2V', retE[:,2,1]                                                ],
                    ['extV', extE[:,2]                                                ],
                    ['Hga', retE[:,3,0]                                                ],
                    ['H2ga', retE[:,3,1]                                                ],
                    ['extga', extE[:,3]                                                ],
                    ['Hgm', retE[:,4,0]                                                ],
                    ['H2gm', retE[:,4,1]                                                ],
                    ['extgm', extE[:,4]                                                ],
                ]
        reactions={}
        for l in setups: # Loop through all species contributions
            reactions[l[0]]=l[1]
        return reactions


    def plot_UE_Erate(self,ratefile,Te,ne,fig=None,ax=0,figsize=(12,13/1.618033),linestyle='-',ylim=(1e-17,1e-11),linewidth=2,labelapp='',savename=None,title=None,pretitle='',figtype='png',plot='loglog',xlim=None,ncol=3,colors=['darkcyan','b','m','r'],eion=5,ediss=10,origUE=False,idx=0,E=0.1,Ton=False,rad=True,Tm=False):
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
        from numpy import log10,array,zeros,sum,ndarray
        from os import mkdir
        from CRUM.reactions import REACTION

        ev=1.602e-19



        rates={}

        try:
            if origUE is True:
                # Read the data from ratefile into dict
                datalist=['H0_depl','H0_create','H2_depl','H2_create','H0_ext','H2_ext']
                self.read_UE(ratefile,rates,datalist=datalist)
                
            else:
                # Read the data from ratefile into dict
                datalist=['Hel','H2el','extel','Hia','H2ia','extia','HV','H2V','extV','Hga','H2ga','extga','Hgm','H2gm','extgm',]
                self.read_UE(ratefile,rates,datalist=datalist)
            # Create custom reaction
            reactions={}
            for r in rates.keys():
                reactions[r]=REACTION(r,'',rates[r],'UE',['',0,'',0,[0,0,0,0]])

        except:
            reactions=self.calc_UE_Erate(Te,ne,E,Ton,rad,Tm)

            '''
            datalist=['Hel','H2el','extel','Hia','H2ia','extia','HV','H2V','extV','Hga','H2ga','extga','Hgm','H2gm','extgm',]
            ret=zeros((len(Te),self.Np,self.Np))
            ext=zeros((len(Te),self.Np))
            retE=zeros((len(Te),5,self.Np))
            extE=zeros((len(Te),5))
            for i in range(len(Te)):
                print(str(i+1)+'/'+str(len(Te)))
                T=Te[i]
                ret[i,:,:],ext[i,:],_=self.crm.gl_crm(*self.crm.M(T,ne,T,ne,E,write=False),Sext=True) # Store to matrix
                U=self.crm.Sgl(T,ne,T,ne,E,rad,Tm,write=False,Ton=Ton) # Don't include erl1/erl2 radiation in the rates as these are handled by UEDGE
                # Include T-losses or not?
                retE[i,:,:]=array([ sum(U[0][0],axis=0), sum(U[1][0],axis=0), sum(U[2][0],axis=0), sum(U[3][0],axis=0), sum(U[4][0],axis=0)])
                #print(len(retE.shape))
                extE[i,:]=array([ sum(U[0][1],axis=0), sum(U[1][1],axis=0), sum(U[2][1],axis=0), sum(U[3][1],axis=0), sum(U[4][1],axis=0) ])


            setups= [   ['Hel', retE[:,0,0]    ],
                        ['H2el', retE[:,0,1]                                                ],
                        ['extel', extE[:,0]                                                ],
                        ['Hia', retE[:,1,0]                                                ],
                        ['H2ia', retE[:,1,1]                                                ],
                        ['extia', extE[:,1]                                                ],
                        ['HV', retE[:,2,0]                                                ],
                        ['H2V', retE[:,2,1]                                                ],
                        ['extV', extE[:,2]                                                ],
                        ['Hga', retE[:,3,0]                                                ],
                        ['H2ga', retE[:,3,1]                                                ],
                        ['extga', extE[:,3]                                                ],
                        ['Hgm', retE[:,4,0]                                                ],
                        ['H2gm', retE[:,4,1]                                                ],
                        ['extgm', extE[:,4]                                                ],
                    ]
            reactions={}
            for l in setups: # Loop through all species contributions
                reactions[l[0]]=l[1]
            '''




        title=['Electron loss','Ion/atom source','Potential source','Radiation source']
        lab=['H0','H2','ext','H2,a','H2,m']
               
        try:
            ax=fig.get_axes()[0]
            oldax=True
        except:
            fig=figure(figsize=figsize)
            oldax=False



        for i in range(4):
            ax=fig.add_subplot(2,2,i+1)

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

            if origUE:

                fac=ev*( ediss*(i==0)+2*eion*(i==1)+(ediss-2*eion)*(i==3)*(ediss>2*eion)-(ediss-2*eion)*(i==2)*(ediss<2*eion))
                # Setup plot type
                if isinstance(Te,list):
                    x=Te
                    y=array([reactions['H2_depl'].rate(T,None,ne=ne)*fac for T in Te])
                    xlabel=(r'Electron temperature [eV]')
                    ex=int(log10(ne))
                    mu=ne/(10**ex) 
                    autotitle=title[i]
                    suptitle=pretitle+r'$\rm{{n_e}}$={}$\times10^{{{}}}$ $\rm{{cm^{{-3}}}}$'.format(mu,ex)
                elif isinstance(ne,list): 
                    x=ne
                    y=array([reactions['H2_depl'].rate(Te,None,ne=n)*fac for n in ne])
                    if idx>3:
                        y=y/ne
                    xlabel=(r'Electron density [$\rm{cm^{-3}}$]')
                    autotitle=title[i]
                    supttitle=pretitle+r' {}, $\rm{{T_e}}$={} [eV]'.format(typ,Te)
                else:
                    print('Reaction rate for process "{}" at Te={} and ne={}'.format(datalist[idx],Te,ne))
                    return reactions[datalist[idx]].rate(Te,None,ne)
                   
                pl(x,((-1)**(i<2))*y,linewidth=linewidth,color='k',label='Orig. UE'+labelapp,linestyle=linestyle)

                ax.set_title(autotitle)
            else:

                if i!=3:
                    ii=3
                else:
                    ii=6

                for j in range(ii):
                    if ii==6:
                        if j>2:
                            jj=3
                        else:
                            jj=j
                    else:
                        jj=j

     
                    # Setup plot type
                    if isinstance(Te,list) or isinstance(Te,ndarray):
                        x=Te
                        try:
                            y=array([reactions[datalist[i*3+j]].rate(T,None,ne=ne)*ev for T in Te])/(ne**(j==2))
                        except:
                            y=reactions[datalist[i*3+j]]*ev
                        xlabel=(r'Electron temperature [eV]')
                        ex=int(log10(ne))
                        mu=ne/(10**ex) 
                        autotitle=title[i]
                        suptitle=pretitle+r'$\rm{{n_e}}$={}$\times10^{{{}}}$ $\rm{{cm^{{-3}}}}$'.format(mu,ex)
                    elif isinstance(ne,list): 
                        x=ne
                        y=array([reactions[datalist[i*3+j]].rate(Te,None,ne=n)*ev for n in ne])
                        if j==2:
                            y=y/ne
                        xlabel=(r'Electron density [$\rm{cm^{-3}}$]')
                        autotitle=title[i]
                        supttitle=pretitle+r' {}, $\rm{{T_e}}$={} [eV]'.format(typ,Te)
                    else:
                        print('Reaction rate for process "{}" at Te={} and ne={}'.format(datalist[idx],Te,ne))
                        return reactions[datalist[idx]].rate(Te,None,ne)

                    if j==0:
                        tot=zeros((len(x),))
                    tot=tot+y

                    pl(x,((-1)**(i==0))*y,color=colors[jj],label=lab[jj]+labelapp,linewidth=linewidth,linestyle=linestyle)

                    ax.set_title(autotitle)

                if (plot not in ['semilogy','loglog']) and (origUE is False):
                    pl(x,((-1)**(i==0))*tot,linewidth=linewidth,color='k',label='Total'+labelapp,linestyle=linestyle)

            if oldax is False: 
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.07,box.width, box.height * 0.95])
            if i==3:
                ax.legend(loc='lower center', bbox_to_anchor=(-0.1, -0.35),frameon=False, ncol=ncol)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(r'$\rm{S_E}$ [$\rm{W/(n_e\cdot V\quad)}$]')
            if xlim is not None:
                ax.set_xlim(xlim)
            ax.set_ylim(ylim)

                    
    
                

        fig.suptitle(suptitle)
        
        if savename is not None:
            try:
                mkdir('output/figs')
            except:
                pass
            fig.savefig('output/figs/{}.{}'.format(savename,figtype),dpi=300,edgecolor=None,format=figtype,bbox_inches='tight')


        fig.show()
        return fig

    def spectrum(self,Te,ne,E=0.1,Ti=None,ni=None,fig=None,units='l',norm=False,figsize=(10,10/1.618033),xlim=None,linewidth=1,split=False,write=False):
        #,Sext=True,False,N=False,fig=None,n0=None,ax=0,figsize=(7+7,7*1.618033),linestyle='-',ylim=None,linewidth=2,labelapp='',qlabel=False,savename=None,title=None,pretitle='',figtype='png',color=None,nuclei=True,Tm=False,rad=True,gl=False,n=None,suptitle=None,ylimE=None,Ton=True):
        ''' Plots atomic and molecular spectra '''
        from matplotlib.pyplot import figure

        try:
            ax=fig.get_axes()[0]
        except:
            fig=figure(figsize=figsize)

        if units=='ev':
            xunit=r'E [eV]'
            if xlim is None:
                xlim=(0,20)
        elif units=='v':
            xunit=r'v [$\mathrm{cm^{-1}}}$]'
            if xlim is None:
                xlim=(1.4e4,1.05e5)
        elif units=='f':
            xunit='Frequency [THz]'
            if xlim is None:
                xlim=(400,3500)
        elif units=='l':
            xunit=r'$\rm{\lambda}$ [nm]'
            if xlim is None:
                xlim=(80,700)
        elif units=='Å':
            xunit=r'$\rm{\lambda}$ [Å]'
            if xlim is None:
                xlim=(800,7000)
            

        if norm is True:
            yunit='Normalized intensity []'
        else:
            yunit=r'Counts [$\rm{s^{-1}}]$'


        data=self.intensity(Te,ne,Ti=None,ni=None,E=E,units=units,norm=norm,write=write)

        if split is True:
            species=['atomic','molecular']
            for i in range(2):
                ax=fig.add_subplot(3,1,i+1)
                for d in range(len(data[i][0,:])):
                    if data[i][0,d]!=0:
                        ax.plot([data[i][0,d],data[i][0,d]],[0,data[i][1,d]],'k',linewidth=linewidth)
                ax.set_ylim(0,1.1*max(data[i][1,:]))
                ax.set_xlim(xlim)
                ax.set_xlabel(xunit)
                ax.set_ylabel(yunit)
                ax.set_title('Total radiation from {} processes'.format(species[i]))


            ax=fig.add_subplot(313)

        else:
            ax=fig.add_subplot(111)

        color=['b','r']
        if norm is True:
            data=self.intensity(Te,ne,Ti=None,ni=None,E=E,units=units,norm=False,write=False)
            n=sum(data[0][1,:])+sum(data[1][1,:])
        else:
            n=1

        for i in range(2):
            ax.plot([],[],'b-')
            ax.plot([],[],'r-')
            for d in range(len(data[i][0,:])):
                if data[i][0,d]!=0:
                    ax.plot([data[i][0,d],data[i][0,d]],[0,data[i][1,d]/n],linewidth=linewidth,color=color[i])
            #print(sum(data[i][0,:]*data[i][1,:]))
            ax.set_ylim((0,1.1*max(max(data[0][1,:]),max(data[1][1,:]))/n))
            ax.set_xlim(xlim)
            ax.set_xlabel(xunit)
            ax.set_ylabel(yunit)
            ax.set_title('Total radiation from molecular processes')
            ax.legend(['Atomic bands','Molceular bands'])

    def plotrate(self,database,name,T,n,E=0.1,logx=True,logy=True,res=200,color='k',linestyle='-',ax=0,figsize=(12,13/1.618033),ylim=[1e-14,1e-6],linewidth=2,savename=None,title='',figtype='png',xlim=None,ncol=3,fig=None):
        from numpy import linspace,logspace,log10
        from matplotlib.pyplot import figure

        try:
            ax=fig.get_axes()[ax]
        except:
            fig=figure(figsize=figsize)
            ax=fig.add_subplot(111)

        if logx is True:
            pf=ax.semilogx

            try:
                x=logspace(log10(T[0]),log10(T[1]))
                y=[self.get_rate(database,name,i,n,E) for i in x]
                xlabel='Temperature [eV]'.format(T)
                titapp=', '*(len(title)>0)+r'n={:.2E} $\rm{{cm^{{-3}}}}$'.format(n)
            except:
                x=logspace(log10(T[0]),log10(T[1]))
                y=[self.get_rate(database,name,T,i,E) for i in x]
                xlabel=r'Density [$\rm{cm^{-3}}$]'
                titapp=', '*(title!='')+r'T={} eV'.format(T)
            if logy is True:
                pf=ax.loglog
        else:
            pf=ax.plot
            try:
                x=linspace(T[0],T[1],res)
                y=[self.get_rate(database,name,i,n,E) for i in x]
                xlabel='Temperature [eV]'.format(T)
                titapp=', '*(len(title)>0)+r'n={:.2E} $\rm{{cm^{{-3}}}}$'.format(n)
            except:
                x=linspace(T[0],T[1],res)
                y=[self.get_rate(database,name,T,i,E) for i in x]
                xlabel=r'Density [$\rm{cm^{-3}}$]'
                titapp=', '*(title!='')+r'T={} eV'.format(T)
            if logy is True:
                pf=ax.semilogy

        pf(x,y,color=color,linestyle=linestyle,linewidth=linewidth,label=database+'.'+name)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$\rm{Rate}$ [$\rm{cm^{3}/s}$]')
        ax.set_title(title+titapp)
        ax.legend(ncol=ncol)

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        

        if savename is not None:
            try:
                mkdir('output/figs')
            except:
                pass
            fig.savefig('output/figs/{}.{}'.format(savename,figtype),dpi=300,edgecolor=None,format=figtype,bbox_inches='tight')


        fig.show()
        return fig









    def steady_state(self,Te,ne,na,nm,diva,divm,srca,srcm,bga,bgm,psor,psorgc,vol,E=0.1,Ti=None,ni=None,Sext=True,gl=False,plot=False):
        ''' Returns and prints a list of radiative lifetimes for each species 
        Te - electron temperature [J]
        ne - electron density [1/m**3]
        srca - total atom source in domain [1/s]
        srcm - total molecule source in domain [1/s]
        bga - background atom source in domain [1/s]
        bgm - background molecule source in domain [1/s]
        psor - ionization current in domain [1/s]
        psorgc - recombination current in domain [1/s]
        vol - volume of domain [m**3]
        '''
        from numpy import log,ones,matmul,polyfit,zeros,exp,sum
        from scipy.optimize import fsolve,minimize
        from scipy.linalg import norm
        from scipy.integrate import solve_ivp
        from matplotlib.pyplot import figure


        def dt(t,n,mat,ext,n0):
            return matmul(mat,n)+ext
        Te/=1.602e-19
        ne/=1e6

    
        # Set volume to be cm**3
        vol*=1e6

        # Set initial density in [1/cm**3]
        n=zeros((len(self.n0()),))
        n[0]=1e-6*na
        n[1]=1e-6*nm

        # Get the rate matrices 
        mat,ext=self.M(Te,ne,Ti,ni,E,write=True) # Get the full rate matrix
        if gl is True:
            mat,ext,n=self.gl_crm(mat,ext,n=n) # Get matrices and eigenvalues/-vectors


        # Set external sources to be [1/cm**3/s]
        ext[0]=(psorgc+diva+srca+bga)/vol
        ext[1]=(divm+srcm+bgm)/vol
        # Normalize ionzation rate to UEDGE values 
        mat[0,0]=(psor/ne)/vol

        
        
        # Simulate to SS
        ss_sol=solve_ivp(lambda x,y: dt(x,y,mat,ext,n),(0,10),n,'LSODA',dense_output=True)

        if plot is True:
            f=figure(figsize=(7,7))
            ax=f.add_subplot(111)
            for i in range(len(n)):
                line='--'
                if 'H2' in self.slist[i]: line='-'

                ax.semilogx(ss_sol.t,ss_sol.y[i,:],line,label=self.slist[i]) 
                if i==0:
                    ax.semilogx(ss_sol.t,ss_sol.y[i,:],'k.',label=self.slist[i]) 
            ax.legend(ncol=3)
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Density [cm**-3]')
        return [ss_sol.y[0,-1]*1e6,ss_sol.y[1,-1]*1e6,1e6*self.totmol(ss_sol.y[:,-1])]
       

















