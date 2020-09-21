# CRUM reaction class reaction.py
# Separated from CRUM.py by holm10
# Changelog
# 200205 - Separated from CRUM.py #holm10
# 200210 - Updated ADAS extrapolation, tidied up code #holm10
 
class REACTION:


    def __init__(self, name, database, coeffs,typ,S,bg,species,Tarr=0):
        ''' Creates an reaction object
            __init__(name,database,reactants,fragments,*keys)
    
            name        -   Reaction ID/name (string)
            database    -   Database where reaction is taken from (string)
            reactants   -   List of strings of reactant ID's/handles
            fragments   -   List of strings of fragment ID's/handles
            coeffs      -   Reaction rate coefficients for interpolation/calculation
            typ         -   Reaction type, either of
                            'UE'    -   UEDGE type interpolation
                            'ADAS'  -   ADAS type interpolation
                            'RATE'  -   EIRENE type interpolation, 1D or 2D depending on shape of coeffs
                            'COEFFICIENT'   -   Spontaneous/radiative transition rate
                            'SIGMA' -   Sawada type cross-section fit

            Optional parameters
            Tarr (0)    -   Temperature array required for ADAS interpolation

            
        '''
        from scipy.interpolate import interp2d,interp1d
        from numpy import array
       
        self.database=database
        self.parseS(*S,bg,species)

        # Store the data required to generate the reaction rates
        self.name=name
        self.coeffs=coeffs
        self.coeffs=array(coeffs)
        self.type=typ
        self.Tarr=array(Tarr)
        self.rec=(not self.e)*(not self.p)
        self.rate=self.pick_rate()

    def parseS(self,r,p,K,bg,species):
        '''Parses Returns list of energies

        Assigns and calculates the energy change of the different comonents
        for a reaction. 

        Parameters
        ----------
        r : 

        '''

        from numpy import ones
        self.r_mult=ones((len(r),))
        self.f_mult=ones((len(p),))
        self.reactants=r
        self.fragments=p
        
        for i in range(len(r)):
            if '*' in r[i]:
                self.r_mult[i]=r[i].split('*')[0]
                self.reactants[i]=r[i].split('*')[1]
        
        for i in range(len(p)):
            if '*' in p[i]:
                self.f_mult[i]=p[i].split('*')[0]
                self.fragments[i]=p[i].split('*')[1]

        try:
            self.K=eval(K) # Convert number to kin. E change
            self.Kdep=False
        except:
            self.K=K
            self.Kdep=True

    
        self.e=('e' in self.reactants)
        self.p=('p' in self.reactants)
        
        self.absorption=False    # Electron absorption
        self.prode=False         # Electron production (not conserving Ee)
        self.radrelax=False      # Radiative relaxation
        self.decay=False         # Non-radiative decay: energy goes into K of prod 

        # Check if reaction is an electron absorption reaction
        if self.e:
            try:
                for x in self.reactants:
                    if x in list(bg): b=x
                if b not in self.fragments:
                    self.absorption=True
                    # Account for two-step processes where electron is absorbed
                    # creating a non-stable particle that decays
                    if len(self.f_mult)>len(self.r_mult): 
                        self.decay=True
            except:
                pass
        
        # Do not consider e-impact ionization impact on e-balance:
        #   all energy supplied to the created electron supplied by
        #   the reactant electron

        else:
            # Proton-impact ionization
            if 'e' in p:
                self.prode=True

        # Radiative relaxation when one particle changes state
        if sum(self.r_mult)==1 and sum(self.f_mult==1): 
            self.radrelax=True
        # Non-radiative decay when one particle becomes many
        if sum(self.r_mult)==1 and sum(self.f_mult)>1:  
            self.decay=True

        # TODO: catches for UEDGE ionization and relaxation

        
        # Total potential of reactants and products
        self.Vr,self.Vp=0,0
        for i in range(len(self.reactants)):
            try:
                self.Vr+=self.r_mult[i]*species[self.reactants[i]]['V']
            except:
                self.Vr+=self.r_mult[i]*bg[self.reactants[i]]['V']
        for i in range(len(self.fragments)):
            try:
                self.Vp+=self.f_mult[i]*species[self.fragments[i]]['V']
            except:
                self.Vp+=self.f_mult[i]*bg[self.fragments[i]]['V']


        
        self.S_r,self.S_V,self.S_g,self.S_e=self.ret0,self.ret0,self.ret0,self.ret0


        # Reactant energy loss
        if self.radrelax is False:
            if self.decay is False:
                if self.Kdep:
                    self.S_r=self.depSr
                else:
                    self.S_r=self.indepSr
            else: # Spontaneous decay
                if self.Kdep:
                    self.S_r=self.depdecaySr
                else:
                    self.S_r=self.indepdecaySr
            self.S_V=self.fSV
        else:
            self.S_V=self.fSV
            self.S_g=self.fSg


        if self.prode is True:
            self.S_e=self.depSe
        else:
            self.S_e=self.ret0


        # Create energy matrix in order to only evaluate conditionals once


        Sl=[self.S_r,self.S_V,self.S_g]
        self.Smat=[]
        
        if self.p:
            if not self.e:
                Sl[0]=self.S_e
                for i in range(4):
                    self.Smat.append([self.ret0,self.ret0])
        
        for i in range(3):
            ext=self.ret0
            if (i!=2):
                self.Smat.append([Sl[i],ext])
            elif (i==2)*(self.database not in ['ADAS','JOHNSON']):
                self.Smat.append([self.ret0,self.ret0])
                self.Smat.append([Sl[i],ext])
            else:
                self.Smat.append([Sl[i],ext])
                self.Smat.append([self.ret0,self.ret0])
                    
        while len(self.Smat)<8:
            for i in range(4):
                self.Smat.append([self.ret0,self.ret0])
            
        
        




 
    def getS(self,r,Te,Ti,Tm,E):
        from numpy import array

        ret=[[column(Te,Ti,Tm,E) for column in row] for row in self.Smat]
        return array(ret)

               

    def ret0(self,*args):
        return 0

    def fSV(self,*args):
        return self.Vp-self.Vr

    def fSg(self,*args):
        return self.Vr-self.Vp

    def depSr(self,Te=0,Ti=0,Tm=0,E=0):
        return self.Vr-self.Vp-eval(self.K)
    
    def indepSr(self,*args):
        return self.Vr-self.Vp-self.K 

    def depdecaySr(self,Te=0,Ti=0,Tm=0,E=0):
        return -eval(self.K)

    def indepdecaySr(self,*args):
        return -self.K

    def depSe(self,Te,*args):
        return Te










    
    def print_reaction(self):
        ''' Returns a string with the reaction 
            print_reaction()
        '''
        ret='{}_{}: '.format(self.database,self.name) # Append reaction ID
        for i in range(len(self.reactants)):    # Loop through all reactants
            ret+=(  self.r_mult[i]!=1)*'{}*'.format(self.r_mult[i])+'''{} 
                '''.format(self.reactants[i])+(i+1!=len(self.reactants))*'+ '

        ret+='=> ' # Add reactants to string
        
        for i in range(len(self.fragments)):    # Loop through all fragmetns
            ret+=(  self.f_mult[i]!=1)*'{}*'.format(self.f_mult[i])+'''{} 
                '''.format(self.fragments[i])+(i+1!=len(self.fragments))*'+ ' 
            # Add fragments to string
        
        return ret 



    def pick_rate(self):
        from scipy.interpolate import interp2d,interp1d
        if self.type=='RATE': 
            return self.EIR_rate

        elif self.type=='COEFFICIENT': 
            return self.coeff_rate

        elif self.type=='SIGMA': 
            return self.SAWADASIGMA_rate

        elif self.type=='ADAS': 
            # Create 1D interpolation function
            self.interpolation=interp1d(self.Tarr,self.coeffs,kind='slinear')   
            return self.ADAS_rate

        elif self.type=='UE': 
            t=[i for i in range(self.coeffs.shape[0])]
            n=[i for i in range(self.coeffs.shape[1])]
            # Create 2D interpolation function
            self.interpolation=interp2d(n,t,self.coeffs,kind='linear') 
            return self.UE_rate

        elif self.type=='APID': 
            x=self.coeffs
            self.apidA=     0.14784*(x==2)      +0.058463*(x==3) 
            self.apidB=[    0.0080871*(x==2)    -0.051272*(x==3), 
                            -0.06227*(x==2)     +0.85310*(x==3), 
                            1.9414*(x==2)       -0.57014*(x==3), 
                            -2.198*(x==2)       +0.76684*(x==3), 
                            0.95894*(x==2),
                            1.133*(x>3),
                            -0.4059*(x>3),
                            0.0714*(x>3)                        ]

            return self.APID_rate
        else: 
            print('Unknown type "{}"'.format(self.type))


    
    def ADAS_rate(self,Te,Ti,E=None,ne=None,omegaj=1):
        T=Te*self.e+Ti*self.p
        if T<self.Tarr[0]:
            Tuse=self.Tarr[0]
            coeff=T/Tuse
        else:
            Tuse=T
            coeff=1
        Tuse=min(Tuse,self.Tarr[-1])
#            T=max(T,self.Tarr[0])
        # TODO: How to deal with extrapolation?
        # TODO: figure out what is implied by the statistical weight omegaj?
        #       set =1 for now
        # Return the rate as calculated from the ADAS fit, per the ADAS manual
        return 2.1716e-8*(1/omegaj)*sqrt(13.6048/Tuse)*self.interpolation(Tuse)


    def EIR_rate(self,Te,Ti,E=None,ne=None,omegaj=1):
        from numpy import log,exp
        T=Te*self.e+Ti*self.p
        ret=0
        if T<0.5:   
            Tuse=0.5
            coeff=T/Tuse
        else:       
            Tuse=T
            coeff=1

        if len(self.coeffs.shape)==2:
            ''' T,E fit '''
            for i in range(9):
                for j in range(9):
                    ret+=self.coeffs[i,j]*(log(Tuse)**i)*(log(E)**j)

        elif len(self.coeffs.shape)==1: 
            ''' T fit '''
            for i in range(9):
                ret+=self.coeffs[i]*(log(Tuse)**i)
                
        else:
            print('''Unknown fit: {}, {}, {}
                    '''.format(self.database,self.name,self.type))

        return coeff*exp(ret)


    def coeff_rate(self,Te,Ti,E=None,ne=None,omegaj=1):
            return self.coeffs
            

    def UE_rate(self,Te,Ti,E=None,ne=None):
        # Turn the temperature and density into log-log variables, 
        # bounded to the limits
        jt=max(0,min(10*(log10(Te+1e-99)+1.2),60))
        jn=max(0,min(2*(log10(ne)-10),15))
        # Interpolate jt
        
        # If the density is being extrapolated, scale accordingly!
        # c=1
        # if self.name in ['RECRAD','IONIZRAD']: c=6.242e11
        return self.interpolation(jn,jt)[0]*\
                    (6.242e11**(self.name in ['RECRAD','IONIZRAD']))

    def SAWADASIGMA_rate(self,Te,Ti,E=None,ne=None):
        from numpy import sqrt,pi,inf,exp
        from scipy.integrate import quad
        T=Te*self.e+Ti*self.p
        # TODO Extend to general species?
        mH2=2*1.6735575e-27 # Assume H2 is reactant 1
        me=9.10938356e-31   # Assume electron is reactant 2
        ev=1.602e-19        # Helper
        VH2=sqrt((2*E*ev)/mH2)  # Get the H2 velocity
        Ve=sqrt((2*T*ev)/me)    # Get the e- velocity
        mr=(me*mH2)/(me+mH2)    # Reduced mass
        vth=sqrt((2*self.coeffs[0]*ev)/mr)  # Thermal CM speed
        
        def sigma(E,Eth,q0,A,Omega,W,gamma,nu): 
                # Calculate sigma from the Sawada fit based on E
            Psi=(nu!=0)*(1-W/E)**nu+(gamma!=0)*(1-(W/E)**gamma) 
                    # Get triplet/singlet Psi
            return (E>=Eth)*q0*(A/W**2)*((W/Eth)**Omega)*Psi    
                    # Perform fit
        
        def R(x,T,coeffs):
            # Integrand function as described in JUEL-3858
            return x*sigma(x*T,*self.coeffs)*exp(-x)

        # Perform integration over velocity space according to JUEL-3858
        return (4/sqrt(pi))*sqrt((T*ev)/(2*me))*\
                    quad(R,0,inf,args=(T,self.coeffs))[0]

    

    def APID_rate(self,Te,Ti,E=None,ne=None):
        from numpy import exp,sqrt,log,pi
        T=Te*self.e+Ti*self.p

        ''' The APID-4 rates are taken from Janev's 1993 IAEA paper, the analytic solutions from Stotler's svlib routine used in DEGAS2 '''
        I=13.6/self.coeffs**2

        def expint(k,p):
            a=[-0.57721566,0.99999193,-0.24991055,0.05519968,
                                                -0.00976004,0.00107857]
            ap=[8.5733287401,18.0590169730,8.6347608925, 0.2677737343 ]
            b=[9.5733223454,25.6329561486,21.0996530827, 3.9584969228 ]

            def en(zn,z):
                return exp(-p)*( 1 + kp/(kp+p)**2 + kp*(kp-2*p)/(kp+p)**4 + 
                            kp*(6*p**2-8*kp*p+kp**2)/(kp+p)**6 )/(kp+p)


            if p<0 or (p==0 and k==1): return None
            elif k*p==0: return (k==0)*exp(-p)/p + (p==0)*(1/(k-1))


            elif p<8 or k==1:
                if p<1:
                    ze1=a[5]
                    for i in [4,3,2,1,0]:
                        ze1=a[i]+ze1*p
                    ze1-=log(p)

                else:
                    znum=1
                    zden=1
                    for i in range(4):
                        znum=ap[i]+znum*p
                        zden=b[i]+zden*p
                    ze1=exp(-p)*znum/(zden*p)

                if k==1:
                    return ze1
                else:
                    zek=ze1
                    for i in range(2,k+1):
                        zek=(exp(-p)-p*zek)/(i-1)
                    return zek
            
            else:
                kp=int(p+0.5)
                if k<kp and k<=10:
                    if kp>10: kp=10
                    zek=en(kp,p)
                    for i in range(kp-1,k-1,-1): zek=(exp(-p) -i*zek)/p
                    return zek
            
                else: return en(k,p)
            return 'Fell througFell throughh'

            #return exp(-p)/(p**k)

        # Use set constants to eval sigma
        if self.coeffs<=3: 
            zarg=I/T

            zeint=[]
            for i in range(1,7): 
                zeint.append(expint(i,zarg))

            zmul=[]
            zmul.append(self.apidB[0] + 2*self.apidB[1] + 3*self.apidB[2] + 
                                        4*self.apidB[3] + 5*self.apidB[4])
            zmul.append(-2*( self.apidB[1] + 3*self.apidB[2] + 6*self.apidB[3]
                                        + 10*self.apidB[4] ))
            zmul.append(3*( self.apidB[2] + 4*self.apidB[3] 
                                        + 10*self.apidB[4] ))
            zmul.append(-4*( self.apidB[3] + 5*self.apidB[4] ))
            zmul.append(5*self.apidB[4])

            zi1=0
            for i in range(5): zi1+=zmul[i]*zeint[i+1]


            zi2=self.apidA*zeint[0]
            zi3=1e-13/(I*T)
            return 6.692e7*sqrt(T) *zi3*(zi1 + zi2)


        # Higher n, evaluate from formulae
        else:

            zn=self.coeffs
            zg0 = 0.9935 + 0.2328/zn - 0.1296/zn**2
            zg1 = -(0.6282 - 0.5598/zn + 0.5299/zn**2) / zn
            zg2 = -(0.3887 - 1.181/zn + 1.470/zn**2) / zn**2
            zmul = 32. / (3. * sqrt(3.) * pi)
            zan = zmul * zn * (zg0/3. + zg1/4. + zg2/5.)  
            zrn = 1.94 * zn**(-1.57)
            zb = (4.0 - 18.63/zn + 36.24/zn**2 - 28.09/zn**3) / zn
            zbn = 2. * zn**2 * (5. + zb) / 3
            zyn = I / T
            zzn = zrn + zyn
            ze1y = expint(1,zyn)
            ze1z = expint(1,zzn)
            zint1 = zan * (ze1y/zyn - ze1z/zzn) 
            zxiy = expint(0,zyn) - 2.*ze1y + expint(2,zyn)
            zxiz = expint(0,zzn) - 2.*ze1z + expint(2,zzn)
            zint2 = (zbn - zan * log(2.*zn**2)) * (zxiy - zxiz)
            zint3 = 1.76e-16 * zn**2 * zyn**2
            ret = 6.692e7 * sqrt(T) * zint3 * (zint1 + zint2)



        return ret
            

