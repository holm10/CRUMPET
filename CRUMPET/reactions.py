# CRUMPET reaction class reaction.py
# Separated from CRUMPET.py by holm10
# Changelog
# 200205 - Separated from CRUMPET.py #holm10
# 200210 - Updated ADAS extrapolation, tidied up code #holm10
 

class Reaction:
    ''' Class containing reaction data

    This class contains data pertaining to a single defined reaction,
    such as energy transfer rates, reactants, products, and reaction
    rates.

    ...
    
    Attributes
    ----------
    name : string
        reaction handle/ID
    database : string
        database where reaction rates are defined
    rate(Te, Ti, ne=None, E=None, omegaj=1) : function
        returns the temperature- and density-dependent reaction rate
        coefficients
        Te : float
            electron temperature in eV
        Ti : float 
            ion temperature in eV
        ne : float, optional (default: None)
            electron temperature in m**-3, only used for type='UE'
        E : float, optional (default: None)
            molecular energy in eV, only used for type='RATE'
        omegaj : float, optional (default: 1)
            statistical weight for ADAS reaction rates
    type : string
        identifier for the kind of reaction being defined:
        'RATE' : EIRENE-type interpolation in one or two dimensions,
                 depending on the shape of coeffs
        'COEFFICIENT' : Spontaneous transition rate, s**-1
        'ADAS' : ADAS-type interpolation
        'SIGMA' : Sawada 1995-type cross-section fit
        'APID' : APID-type ionization rate
        'UE' : UEDGE-type interpolation
    reactants : list of strings
        a list of the reactant species handles
    fragments : list of strings
        a list of the fragment species handles
    Vp : float
        the total potential of the products
    Vr : float
        the total potential of the reactants
    K : float/string
        the kinetic energy transferred to the procucts in the reaction:
        if a species is born in the reaction, they are assigned the 
        local temperature of that species, which is evaluated upon
        calls to getS
    Kdep : boolean
        marker determining wheter K is dependent on the species 
        temperature, i.e. whether a new species is born
    S_V(*args) : function
        Returns the change in potential for the reaction calculated
        from the potential of the reactions and products
    S_e(Te,*args) : function
        Returns the electron energy increase for the reaction. Only
        relevant for ionization reactions. Te is the electron 
        temperature in eV
    S_g(*args) : function
        Returns the gamma energy for the reaction calculated from the
        potential difference between the reactants and products. Only
        nonzero if the reaction is radiative decay
    S_r(Te,Ti,Tm) : function
        Returns the background reactant energy loss, calculated from
        the net energy balance
    Smat : 2D list
        Smat has len=8 along axis=0 and len=2 along axis=1. The first 4
        entries along axis=0 are: reactant energy loss, potential energy
        change, radiated energy in atomic lines, and radiated energy in 
        molecular lines, all for electron-impact processes. The next 4 
        entries are the corresponding energy terms but for proton-impact
        processes. The first entries along axis=1 are the internal terms
        and the second due to external sinks and sources (currently 0)
    Tarr : ndarray
        A vector containing the electron temperature points of ADAS fits
        only relevant for ADAS rates, otherwise not used
    absorption : boolean
        marker determining if the reaction is an electron-absorption
        reaction
    coeffs : ndarray
        the reaction rate coefficients. The shape and values depend on
        the type of reaction. Used by rate() to calculate the reaction
        rates
    decay : boolean
        marker whether the reaction is non-radiative decay, e.g. decay
        from an unstable state into dissociation products
    radrelax : boolean
        marker whether the reaction is radiative relaxation
    prode : boolean
        marker whether an electron is created in the reaction, i.e. the
        reaction does not conserve electron energy
    e : boolean
        marker whether the reaction is an electron-impact reaction or
        not
    p : boolean 
        marker whether the reaction is and proton-impact reaction or not
    f_mult : int
        list of fragment mutipliers, same length as fragments
    r_mult : int
        list of product mutipliers, same length as products
    
    Methods
    ------- 
    '''


    def __init__(self, database, rtype, coeffs, typ, data, bg, species,
            isotope='H', mass=1, Tarr=0):
        ''' 
        Parameters 
        ----------
        name : string
            reaction handle/ID
        database : string
            database where reaction rates are defined
        rtype : string
            database fit type according to EIRENE definitions
        coeffs : ndarray
            the reaction rate coefficients. The shape and values depend
            on the type of reaction. Used by rate() to calculate the
            reaction rates
        type : string
            identifier for the kind of reaction being defined:
            'RATE' : EIRENE-type interpolation in one or two dimensions,
                     depending on the shape of coeffs
            'COEFFICIENT' : Spontaneous transition rate, s**-1
            'ADAS' : ADAS-type interpolation
            'SIGMA' : Sawada 1995-type cross-section fit
            'APID' : APID-type ionization rate
            'UE' : UEDGE-type interpolation
        S : dict
            dict containing the reaction information, requiring the 
            following keys:
            'reactants' : list of strings
                reactant handles
            'fragments' : list of strings
                fragment handles
            'K' : string
                string of kinetic energy information
        bg : dict
            Crm.bg, dict of background species and potentials
        species : dict
            Crm.species, dict of CRM species and potentials
        Tarr : list, optional (default: 0)
            Array of ADAS temperature point. Only used when typ=ADAS
        isotope : string
            The handle for the isotope being used. Used to identify 
            atomic and molecular species
        mass : float
            The isotope mass in AMU, used for calculations in Crm 
        '''
        from scipy.interpolate import interp2d, interp1d
        from numpy import array

        self.parseS(data, bg, species, database, isotope)
        self.isotope = isotope
        self.mass = mass
        self.coeffs = coeffs
        self.coeffs = array(coeffs)
        self.type = typ
        self.fittype = rtype
        self.Tarr = array(Tarr)
        self.rate = self.pick_rate()


    def parseS(self, data, bg, species, database, isotope):
        ''' Derives the reaction information
        
        Assigns the majority of the object attributes, including
        the energy sinks/sources for the different components       

        Parameters
        ----------
        r : list of strings
            reactant handles
        p : list of strings
            fragment handles
        K : string
            string of kinetic energy information
        bg : dict
            Crm.bg, dict of background species and potentials
        species : dict
            Crm.species, dict of CRM species and potentials

        Returns
        -------
        None
        '''
        from numpy import ones

        r = data['reactants']
        p = data['fragments']
        K = data['K']
        self.r_mult = ones((len(r),))
        self.f_mult = ones((len(p),))
        self.reactants = r
        self.fragments = p
        for i in range(len(r)):
            if '*' in r[i]:
                self.r_mult[i] = r[i].split('*')[0]
                self.reactants[i] = r[i].split('*')[1]
        
        for i in range(len(p)):
            if '*' in p[i]:
                self.f_mult[i] = p[i].split('*')[0]
                self.fragments[i] = p[i].split('*')[1]
        try:
            self.K = eval(K) # Convert number to kin. E change
            self.Kdep = False
        except:
            self.K = K
            self.Kdep = True
        self.e = ('e' in self.reactants)
        self.p = ('p' in self.reactants)
        self.absorption = False    # Electron absorption
        self.prode = False         # Electron production (not conserving Ee)
        self.radrelax = False      # Radiative relaxation
        self.decay = False         # Non-radiative decay: energy goes into 
                                   # K of prod 
        # TODO: how to treat ionization reactions? Included in K? Or not?
        # Check if reaction is an electron absorption reaction
        if self.e:
            try:
                for x in self.reactants:
                    if x in list(bg): 
                        b = x
                if b not in self.fragments:
                    self.absorption = True
                    # Account for two-step processes where electron is absorbed
                    # creating a non-stable particle that decays
                    if len(self.f_mult) > len(self.r_mult): 
                        self.decay = True
            except:
                pass
        # Do not consider e-impact ionization impact on e-balance:
        #   all energy supplied to the created electron supplied by
        #   the reactant electron
        else:
            # Proton-impact ionization
            if 'e' in p:
                self.prode = True
        # Radiative relaxation when one particle changes state
        if sum(self.r_mult) == 1 and sum(self.f_mult == 1): 
            self.radrelax=True
        # Non-radiative decay when one particle becomes many
        if sum(self.r_mult) == 1 and sum(self.f_mult) > 1:  
            self.decay = True
        # TODO: catches for UEDGE ionization and relaxation
        # Total potential of reactants and products
        self.Vr = 0
        self.Vp = 0
        for i in range(len(self.reactants)):
            try:
                self.Vr += self.r_mult[i]*species[self.reactants[i]]['V']
            except:
                self.Vr += self.r_mult[i]*bg[self.reactants[i]]['V']
        for i in range(len(self.fragments)):
            try:
                self.Vp += self.f_mult[i]*species[self.fragments[i]]['V']
            except:
                self.Vp += self.f_mult[i]*bg[self.fragments[i]]['V']
        self.S_r = self.ret0
        self.S_V = self.ret0
        self.S_g = self.ret0
        self.S_e = self.ret0
        # Reactant energy loss
        if self.radrelax is False:
            if self.decay is False:
                if self.Kdep:
                    self.S_r = self.depSr
                else:
                    self.S_r = self.indepSr
            else: # Spontaneous decay
                if self.Kdep:
                    self.S_r = self.depdecaySr
                else:
                    self.S_r = self.indepdecaySr
            self.S_V = self.fSV
        else:
            self.S_V = self.fSV
            self.S_g = self.fSg
        if self.prode is True:
            self.S_e = self.depSe
        else:
            self.S_e = self.ret0
        # Create energy matrix in order to only evaluate conditionals once

        # TODO: Verify external energy change - when is it relevant
        Sl = [self.S_r, self.S_V, self.S_g]
        self.Smat = []
        if self.p:
            if not self.e:
                Sl[0] = self.S_e
                for i in range(4):
                    self.Smat.append([self.ret0, self.ret0])
        for i in range(3):
            ext = self.ret0
            if (i != 2):
                self.Smat.append([Sl[i], ext])
            elif (i == 2)*('{}2'.format(isotope) in ''.join(self.reactants)):
                self.Smat.append([self.ret0, self.ret0])
                self.Smat.append([Sl[i], ext])
            else:
                self.Smat.append([Sl[i], ext])
                self.Smat.append([self.ret0, self.ret0])
        while len(self.Smat) < 8:
            for i in range(4):
                self.Smat.append([self.ret0, self.ret0])
            

    def getS(self, Te, Ti, Tm, E):
        ''' Evaluates Smat at given Te, Ti, Tm, and E
        
        Parameters
        ----------
        Te : float
            electron temperature in eV
        Ti : float
            ion temperature in eV
        Tm : float
            molecular temperature in eV
        E : float
            assumed molecular rest energy in eV

        Returns
        -------
        ndarray
            Evaluated ndarray of Smat
        '''
        from numpy import array,sum
        if Tm == 0:
            Tm=E
        S=array([[column(Te, Ti, Tm) for column in row] for row in self.Smat])
        ret = array([
                [S[0] + S[4]],
                [-sum(S, axis=0)],
                [S[1] + S[5]],
                [S[2] + S[6]],
                [S[3] + S[7]] ],
                )[:,0,:]
        return ret
               

    def ret0(self, *args):
        ''' Zero-function'''
        return 0


    def fSV(self, *args):
        ''' Potential difference '''
        return self.Vp - self.Vr


    def fSg(self, *args):
        ''' Gamma energy'''
        return self.Vr - self.Vp


    def depSr(self, Te=0, Ti=0, Tm=0):
        ''' Temperature-dependent reactant sink '''
        return self.Vr - self.Vp - eval(self.K)

    
    def indepSr(self, *args):
        ''' Temperautre-independent reactant sink '''
        return self.Vr - self.Vp - self.K


    def depdecaySr(self, Te=0, Ti=0, Tm=0):
        ''' Temperature-dependent fragmentation '''
        return -eval(self.K)


    def indepdecaySr(self, *args):
        ''' Temperature-independent fragmentation '''
        return self.K


    def depSe(self, Te, *args):
        ''' Electron energy gain '''
        return Te

    
    def print_reaction(self, database, name):
        ''' Returns formatted a string with the reaction '''
        ret = '{}_{}_{}: '.format(database, self.fittype, name) # Append reaction ID
        for i in range(len(self.reactants)):    # Loop through all reactants
            ret += (self.r_mult[i] != 1)*'{}*'.format(self.r_mult[i])\
                    + '{}'.format(self.reactants[i])\
                    + (i+1 != len(self.reactants))*'+ '
        ret += ' => ' # Add reactants to string
        for i in range(len(self.fragments)):    # Loop through all fragmetns
            ret += (self.f_mult[i] != 1)*'{}*'.format(self.f_mult[i])\
                    + '{}'.format(self.fragments[i])\
                    + (i+1 != len(self.fragments))*'+ ' 
            # Add fragments to string
        
        return ret 


    def pick_rate(self):
        ''' Initialization of the rate attribute '''
        # TODO: Dynamically choose between ni and ne in reaction
        from scipy.interpolate import interp2d,interp1d
        if self.type == 'RATE': 
            return self.EIR_rate
        elif self.type == 'COEFFICIENT': 
            return self.coeff_rate
        elif self.type == 'SIGMA': 
            return self.SAWADASIGMA_rate
        elif self.type == 'ADAS': 
            # Create 1D interpolation function
            self.interpolation = interp1d(
                    self.Tarr, self.coeffs, kind='slinear')   
            return self.ADAS_rate
        elif self.type == 'UE': 
            t = [i for i in range(self.coeffs.shape[0])]
            n = [i for i in range(self.coeffs.shape[1])]
            # Create 2D interpolation function
            self.interpolation = interp2d(n, t, self.coeffs, kind='linear') 
            return self.UE_rate
        elif self.type == 'APID': 
            x = self.coeffs
            self.apidA =  0.14784*(x == 2) + 0.058463*(x == 3) 
            self.apidB = [
                            0.0080871*(x == 2) - 0.051272*(x == 3), 
                            -0.06227*(x == 2) + 0.85310*(x == 3), 
                            1.9414*(x == 2) - 0.57014*(x == 3), 
                            -2.198*(x == 2) + 0.76684*(x == 3), 
                            0.95894*(x == 2),
                            1.133*(x > 3),
                            -0.4059*(x > 3),
                            0.0714*(x > 3)
                        ]
            return self.APID_rate
        else: 
            print('Unknown type "{}"'.format(self.type))

    
    def ADAS_rate(self, Te, Ti, omegaj=1, **kwargs):
        ''' Function returning the ADAS rate for T '''
        T = Te*self.e + Ti*self.p
        if T < self.Tarr[0]:
            Tuse = self.Tarr[0]
            coeff = T/Tuse
        else:
            Tuse = T
            coeff = 1
        Tuse = min(Tuse, self.Tarr[-1])
#            T=max(T,self.Tarr[0])
        # TODO: How to deal with extrapolation?
        # TODO: figure out what is implied by the statistical weight omegaj?
        #       set =1 for now
        # Return the rate as calculated from the ADAS fit, per the ADAS manual
        return 2.1716e-8*(1/omegaj)*sqrt(13.6048/Tuse)*self.interpolation(Tuse)


    def EIR_rate(self, Te, Ti, E=None, ne=None, **kwargs):
        ''' Function returning the EIRENE rate for T '''
        from numpy import log, exp
        T = Te*self.e + Ti*self.p
        ret = 0
        if T < 0.5:   
            Tuse = 0.5
            coeff = T/Tuse
        else:       
            Tuse = T
            coeff = 1
        if self.fittype == 'H.0':
            print('Potential: TBD')
            return
        elif self.fittype == 'H.1':
            print('Cross-section vs energy: TBD')
            return
        elif self.fittype == 'H.2':
            # Rate coefficient vs temperature
            for i in range(9):
                ret += self.coeffs[i]*(log(Tuse)**i)    
        elif self.fittype == 'H.3':
            # Rate coefficient vs temperature and energy
            for i in range(9):
                for j in range(9):
                    ret += self.coeffs[i,j]*(log(Tuse)**i)*(log(E)**j)
        elif self.fittype == 'H.4':
            # Rate coefficient vs temperature and density
            for i in range(9):
                for j in range(9):
                    ret += self.coeffs[i,j]*(log(Tuse)**i)*(log(ne*1e-8)**j)
        elif self.fittype == 'H.5':
            print('Momentum-weighted rates vs. temperature, not in use')
            return
        elif self.fittype == 'H.6':
            print('Momentum-weighted rates vs. temperaturei and energy: TBD')
            return
        elif self.fittype == 'H.7':
            print('Momentum-weighted rates vs. temperature and density, not in use')
            return
        elif self.fittype == 'H.8':
            print('Energy-weighted rates vs. temperature: TBD')
            return
        elif self.fittype == 'H.9':
            print('Energy-weighted rates vs. temperature and energy, not in use')
            return
        elif self.fittype == 'H.10':
            print('Energy-weighted rates vs. temperature and density: TBD')
            return
        elif self.fittype == 'H.11':
            print('Other data: TBD')
            return
        elif self.fittype == 'H.12':
            print('Other data: TBD')
            return
        else:
            print('Unknown fit: {}, {}, {}'.format(
                    self.database,self.name,self.type))
        return coeff*exp(ret)


    def coeff_rate(self, *args, **kwargs):
        return self.coeffs
            

    def UE_rate(self, Te, *args, E=None, ne=None, **kwargs):
        ''' Function returning the UEDGE rate for T and n '''
        from numpy import log10
        # Turn the temperature and density into log-log variables, 
        # bounded to the limits
        jt = max(0,min(10*(log10(Te + 1e-99)+1.2), 60))
        jn = max(0,min(2*(log10(ne) - 10), 15))
        # Interpolate jt
        
        # If the density is being extrapolated, scale accordingly!
        # c=1
        # if self.name in ['RECRAD','IONIZRAD']: c=6.242e11
        return self.interpolation(jn, jt)[0]#*\
        #        (6.242e11**(self.name in ['RECRAD', 'IONIZRAD']))


    def SAWADASIGMA_rate(self, Te, Ti, **kwargs):
        ''' Function returning the Sawada-Sigma rate for T '''
        from numpy import sqrt,pi,inf,exp
        from scipy.integrate import quad
        T = Te*self.e + Ti*self.p
        # TODO Extend to general species?
        mm = self.mass*2*1.6735575e-27 
        me = 9.10938356e-31   # Assume electron is reactant 2
        ev = 1.602e-19        # Helper
        #VH2=sqrt((2*E*ev)/mm)  # Get the H2 velocity
        Ve = sqrt((2*T*ev)/me)    # Get the e- velocity
        mr = (me*mm)/(me+mm)    # Reduced mass
        vth = sqrt((2*self.coeffs[0]*ev)/mr)  # Thermal CM speed
        
        def sigma(E, Eth, q0, A, Omega, W, gamma, nu): 
                # Calculate sigma from the Sawada fit based on E
            Psi = (nu != 0)*(1-W/E)**nu+(gamma != 0)*(1 - (W/E)**gamma) 
                    # Get triplet/singlet Psi
            return (E >= Eth)*q0*(A/W**2)*((W/Eth)**Omega)*Psi    
                    # Perform fit
        
        def R(x, T, coeffs):
            # Integrand function as described in JUEL-3858
            return x*sigma(x*T, *self.coeffs)*exp(-x)

        # Perform integration over velocity space according to JUEL-3858
        return (4/sqrt(pi))*sqrt((T*ev)/(2*me))*\
                quad(R, 0, inf, args=(T, self.coeffs))[0]
    

    def APID_rate(self, Te, Ti, **kwargs):
        ''' Function returning the APID rate for T '''
        from numpy import exp, sqrt, log, pi
        T = Te*self.e + Ti*self.p
        ''' The APID-4 rates are taken from Janev's 1993 IAEA paper, the 
            analytic solutions from Stotler's svlib routine used in DEGAS2 '''
        I = 13.6/self.coeffs**2

        def expint(k, p):
            a = [-0.57721566, 0.99999193, -0.24991055, 0.05519968, -0.00976004,
                    0.00107857]
            ap = [8.5733287401, 18.0590169730, 8.6347608925, 0.2677737343]
            b = [9.5733223454, 25.6329561486, 21.0996530827, 3.9584969228]

            def en(zn, z):
                return exp(-p)*( 1 + kp/(kp+p)**2 + kp*(kp - 2*p)/(kp + p)**4 +
                        kp*(6*p**2 - 8*kp*p+kp**2)/(kp + p)**6 )/(kp + p)

            if p < 0 or (p == 0 and k == 1): 
                return None
            elif k*p == 0: 
                return (k == 0)*exp(-p)/p + (p == 0)*(1/(k-1))

            elif p < 8 or k == 1:
                if p < 1:
                    ze1 = a[5]
                    for i in [4,3,2,1,0]:
                        ze1 = a[i] + ze1*p
                    ze1 -= log(p)

                else:
                    znum = 1
                    zden = 1
                    for i in range(4):
                        znum = ap[i] + znum*p
                        zden = b[i] + zden*p
                    ze1 = exp(-p)*znum/(zden*p)

                if k == 1:
                    return ze1
                else:
                    zek = ze1
                    for i in range(2, k+1):
                        zek = (exp(-p) - p*zek)/(i - 1)
                    return zek
            
            else:
                kp = int(p + 0.5)
                if k < kp and k <= 10:
                    if kp > 10: 
                        kp = 10
                    zek = en(kp, p)
                    for i in range(kp - 1, k - 1, -1): 
                        zek=(exp(-p) - i*zek)/p
                    return zek
            
                else: return en(k, p)
            return 'Fell through'

        # Use set constants to eval sigma
        if self.coeffs <= 3: 
            zarg = I/T
            zeint = []
            for i in range(1, 7): 
                zeint.append(expint(i, zarg))
            zmul = []
            zmul.append(self.apidB[0] + 2*self.apidB[1] + 3*self.apidB[2] + 
                    4*self.apidB[3] + 5*self.apidB[4])
            zmul.append(-2*(self.apidB[1] + 3*self.apidB[2] + 6*self.apidB[3]
                    + 10*self.apidB[4] ))
            zmul.append(3*(self.apidB[2] + 4*self.apidB[3] 
                    + 10*self.apidB[4] ))
            zmul.append(-4*(self.apidB[3] + 5*self.apidB[4] ))
            zmul.append(5*self.apidB[4])
            zi1 = 0
            for i in range(5):
                zi1 += zmul[i]*zeint[i + 1]
            zi2 = self.apidA*zeint[0]
            zi3 = 1e-13/(I*T)
            return 6.692e7*sqrt(T)*zi3*(zi1 + zi2)

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
            ze1y = expint(1, zyn)
            ze1z = expint(1, zzn)
            zint1 = zan * (ze1y/zyn - ze1z/zzn) 
            zxiy = expint(0, zyn) - 2.*ze1y + expint(2, zyn)
            zxiz = expint(0, zzn) - 2.*ze1z + expint(2, zzn)
            zint2 = (zbn - zan * log(2.*zn**2)) * (zxiy - zxiz)
            zint3 = 1.76e-16 * zn**2 * zyn**2
            ret = 6.692e7 * sqrt(T) * zint3 * (zint1 + zint2)

        return ret
            

