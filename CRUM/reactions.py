# CRUM reaction class reaction.py
# Separated from CRUM.py by holm10
# Changelog
# 200205 - Separated from CRUM.py #holm10
# 200210 - Updated ADAS extrapolation, tidied up code #holm10
 
class REACTION:


    def __init__(self, name, database, reactants, fragments, coeffs,typ,S,Tarr=0):
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
        from numpy import ones,array
        from scipy.interpolate import interp2d,interp1d
        
        # Store the data required to generate the reaction rates
        self.name=name
        self.coeffs=coeffs
        self.database=database
        self.reactants=[r.strip() for r in reactants]
        self.fragments=[f.strip() for f in fragments]
        self.coeffs=array(coeffs)
        self.type=typ
        self.Tarr=array(Tarr)
        for i in range(4):
            if isinstance(S[i],str): S[i]=S[i][4:]

        self.S_r=S[0]
        self.S_g=S[1]
        self.S_V=S[2]
        self.S_e=S[3]

        # Setup multilplier arrays
        self.r_mult=ones((len(reactants),))
        self.f_mult=ones((len(fragments),))
        
        # Separate coefficients from reactant species
        for i in range(len(self.reactants)):    # Loop through reactants
            if '*' in self.reactants[i]:    # If there is a multiplier, separate
                self.f_mult[i]=float(reactants[i].split('*')[0])    # Update multiplier
                self.reactants[i]=reactants[i].split('*')[1].strip()    # Remove mult from reactant string

        # Separate coefficients from fragment species
        for i in range(len(self.fragments)):    # Loop through fragments
            if '*' in self.fragments[i]:    # If there is a multiplier, separate
                self.f_mult[i]=float(fragments[i].split('*')[0])    # Update the multiplier
                self.fragments[i]=fragments[i].split('*')[1].strip()    # Remove the multiplier from fragment string


        ''' Make interpolation object if necessary '''
        if self.type=='UE': # UEDGE interpolator
            t=[i for i in range(self.coeffs.shape[0])]
            n=[i for i in range(self.coeffs.shape[1])]

            self.interpolation=interp2d(n,t,self.coeffs) # Create 2D interpolation function

        if self.type=='ADAS': # ADAS interpolator
            self.interpolation=interp1d(self.Tarr,self.coeffs,kind='slinear')   # Create 1D interpolation function


    
    def print_reaction(self):
        ''' Returns a string with the reaction 
            print_reaction()
        '''
        ret='' # Initialize string
        ret+='{}_{}: '.format(self.database,self.name) # Append reaction ID
        for i in range(len(self.reactants)):    # Loop through all reactants
            ret+=(self.r_mult[i]!=1)*'{}*'.format(self.r_mult[i])+'{} '.format(self.reactants[i])+(i+1!=len(self.reactants))*'+ '
        ret+='=> ' # Add reactants to string
        for i in range(len(self.fragments)):    # Loop through all fragmetns
            ret+=(self.f_mult[i]!=1)*'{}*'.format(self.f_mult[i])+'{} '.format(self.fragments[i])+(i+1!=len(self.fragments))*'+ ' # Add fragments to string
        
        return ret 

            
    def rate(self,Te,Ti,E=None,ne=None,omegaj=1):
        ''' Returns the reaction rate at the specified temperature
            rate(Te,Ti,*keys)

            Te          -   Electron temperature for evaluation, used if electron reaction [eV]
            Ti          -   Ion temperature for evaluation, used if proton impact [eV]
            
            Optional parameters
            E (None)    -   Target particle energy, used if proton impact [eV]
            ne (None)   -   Electron density, used for UEDGE rates [cm**3]
            omegaj (1)  -   Statistical weight of ADAS rates 
        '''
        from numpy import log,exp,sqrt,pi,inf,array,log10
        from scipy.integrate import quad

        # Find reactant species
        if 'e' in self.reactants:   T=Te # Electron-mediated reaction, use Te
        elif 'p' in self.reactants: T=Ti # Proton impact reaction, use Ti
        else: T=0   # Photoemission, T not used

        # Get rate based on self.type
        if self.type=='RATE':
            ''' We have an EIRENE polynomial fit '''
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
            ''' SAWADA cross-section '''
            # TODO Extend to general species?
            mH2=2*1.6735575e-27 # Assume H2 is reactant 1
            me=9.10938356e-31   # Assume electron is reactant 2
            ev=1.602e-19        # Helper
            VH2=sqrt((2*E*ev)/mH2)  # Get the H2 velocity
            Ve=sqrt((2*T*ev)/me)    # Get the e- velocity
            mr=(me*mH2)/(me+mH2)    # Reduced mass
            vth=sqrt((2*self.coeffs[0]*ev)/mr)  # Thermal CM speed
            
            def sigma(E,Eth,q0,A,Omega,W,gamma,nu): # Calculate sigma from the Sawada fit based on E
                Psi=(nu!=0)*(1-W/E)**nu+(gamma!=0)*(1-(W/E)**gamma) # Get triplet/singlet Psi
                return (E>=Eth)*q0*(A/W**2)*((W/Eth)**Omega)*Psi    # Perform fit
            
            def R(x,T,coeffs):
                # Integrand function as described in JUEL-3858
                return x*sigma(x*T,*self.coeffs)*exp(-x)

            # Perform integration over velocity space according to JUEL-3858
            return (4/sqrt(pi))*sqrt((T*ev)/(2*me))*quad(R,0,inf,args=(T,self.coeffs))[0]


        elif self.type=='ADAS':
            ''' ADAS fit '''
            T=max(T,self.Tarr[-1])
            # TODO: How to deal with extrapolation?
            # TODO: figure out what is implied by the statistical weight omegaj - set =1 for now
            # Return the rate as calculated from the ADAS fit, per the ADAS manual
            return 2.1716e-8*(1/omegaj)*sqrt(13.6048/T)*self.interpolation(T) 
        
        elif self.type=='UE':
            ''' UEDGE fit '''
            # Turn the temperature and density into log-log variables, bounded to the limits
            jt=max(0,min(10*(log10(Te)+1.2),60))
            jn=max(0,min(2*(log10(ne)-10),15))
            # Interpolate jt
            
            c=1
            if self.name in ['RECRAD','IONIZRAD']: c=6.242e11
            return self.interpolation(jn,jt)[0]*c
                    
        else:
            print('Unknown type "{}"'.format(self.type))
    



