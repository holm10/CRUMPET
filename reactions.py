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
    products : list of strings
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
    p_mult : int
        list of fragment mutipliers, same length as products
    e_mult : int
        list of product mutipliers, same length as products
    
    Methods
    ------- 
    '''


    def __init__(self, database, rtype, name, data, coeffs=None, bg=None,
                    species=None, isotope='H', mass=1, Tarr=0, scale=1,
                    sourcefile=None):
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
            'products' : list of strings
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
        from numpy import array, log10

        self.database = database.upper()
        self.isotope = isotope.upper()
        self.name = name
        self.type = rtype.upper()
        self.scale = scale
        self.mass = mass
        self.sourcefile = sourcefile
        self.K = '0'
        self.tag = '{} {} {}'.format(self.database, self.type, self.name)

        reaction = data.pop(0)
        self.educts = [x.strip() for x in \
                            reaction.split(' > ')[0].split(' + ')]
        self.products = [x.strip() for x in \
                            reaction.split(' > ')[1].split(' + ')]
        self.e_mult, self.p_mult = [],[]

        for i in range(len(self.educts)):
            if '*' in self.educts[i]:
                m, s = self.educts[i].split('*')
                self.educts[i] = s.strip()
                self.e_mult.append(float(m))
            else:
                self.e_mult.append(1)

        for i in range(len(self.products)):
            if '*' in self.products[i]:
                m, s = self.products[i].split('*')
                self.products[i] = s.strip()
                self.p_mult.append(float(m))
            else:
                self.p_mult.append(1)

        self.e_mult = array(self.e_mult)
        self.p_mult = array(self.p_mult)
       # print(database, rtype, name, data, coeffs)
        # Make indices array if the reaction is a direct EIRENE rate taken 
        # from a database
        # Read user input rates
        self.coeffs = coeffs
        if self.database == 'USER':
            self.coeffs = []
            if self.type == 'H.2':
                self.coeffs = list(map(float, data.pop(0).split()))
            elif self.type in ['H.3', 'H.4']:
                for i in range(9):
                    self.coeffs.append(list(map(int, data.pop(0).split())))
                self.coeffs = array(self.coeffs)
            elif self.type.upper() == 'RELAXATION':
                self.coeffs = float(data.pop(0))
            elif self.type.upper() == 'COEFFICIENT':
                self.coeffs = float(data[0])
            elif self.type.upper() == 'SIGMA':
                self.coeffs = [float(x) for x in data.pop(0).split()]
            elif self.type.upper() == 'INTERPOLATION':
                pass
            else:
                print('Reaction database "{}" and type "{}"'
                        ' not recognized!'.format(self.database, self.type))

        if self.type.upper() == 'INTERPOLATION':
            # Read the number of data points defined
            self.coeffs = []
            for i in range(int(data.pop(0))):
                self.coeffs.append([float(x) for x in data.pop(0).strip().split()])
            self.coeffs=log10(array(self.coeffs))

        # If user-defined, get it from data
#        if database.upper() == 'USER':
#            self.coeffs = data.pop(0) 

        for s in data:
            if ('K=' in s) or ('K =' in s):
                self.K = s.split('=')[1].strip()
#                try:
#                    self.K = eval(self.K)
#                except:
#                    pass





        self.Tarr = array(Tarr)
        self.parseS(bg, species, isotope)
        self.rate = self.pick_rate()

        '''
        print('==============')
        print(self.database, self.type, self.name)
        print(self.coeffs)
        print(self.rate(1,1,ne=1e12,E=0.1))
        print(self.rate(10,10,ne=1e13,E=1))'''


    def parseS(self, bg, species, isotope):
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

        self.e = ('e' in self.educts)
        self.p = ('p' in self.educts)
        self.absorption = False    # Electron absorption
        self.prode = False         # Electron production (not conserving Ee)
        self.radrelax = False      # Radiative relaxation
        self.decay = False         # Non-radiative decay: energy goes into 
                                   # K of prod 
        # TODO: how to treat ionization reactions? Included in K? Or not?
        # Check if reaction is an electron absorption reaction
        if self.e:
            try:
                for x in self.educts:
                    if x in list(bg): 
                        b = x
                if b not in self.products:
                    self.absorption = True
                    # Account for two-step processes where electron is absorbed
                    # creating a non-stable particle that decays
                    if len(self.p_mult) > len(self.e_mult): 
                        self.decay = True
            except:
                pass
        # Do not consider e-impact ionization impact on e-balance:
        #   all energy supplied to the created electron supplied by
        #   the reactant electron
        else:
            # Proton-impact ionization
            if 'e' in self.products:
                self.prode = True
        # Radiative relaxation when one particle changes state
        if sum(self.e_mult) == 1 and sum(self.p_mult == 1): 
            self.radrelax=True
        # Non-radiative decay when one particle becomes many
        if sum(self.e_mult) == 1 and sum(self.p_mult) > 1:  
            self.decay = True
        # TODO: catches for UEDGE ionization and relaxation
        # Total potential of reactants and products
        self.Vr = 0
        self.Vp = 0
        for i in range(len(self.educts)):
            try:
                self.Vr += self.e_mult[i]*species[self.educts[i]]['V']
            except:
                self.Vr += self.e_mult[i]*bg[self.educts[i]]['V']
        for i in range(len(self.products)):
            try:
                self.Vp += self.p_mult[i]*species[self.products[i]]['V']
            except:
                self.Vp += self.p_mult[i]*bg[self.products[i]]['V']
        self.S_r = self.ret0
        self.S_V = self.ret0
        self.S_g = self.ret0
        self.S_e = self.ret0
        # Reactant energy loss
        if self.radrelax is False:
            if self.decay is False:
                self.S_r = self.Sr
            else: # Spontaneous decay
                self.S_r = self.decaySr
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
            elif (i == 2)*('{}2'.format(isotope) in ''.join(self.educts)):
                self.Smat.append([self.ret0, self.ret0])
                self.Smat.append([Sl[i], ext])
            else:
                self.Smat.append([Sl[i], ext])
                self.Smat.append([self.ret0, self.ret0])
        while len(self.Smat) < 8:
            for i in range(4):
                self.Smat.append([self.ret0, self.ret0])
            

    def getS(self, Te, Ti, Tm, E, marker=False):
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
        if marker is True:
            ret[ret<0]=0
            ret[ret!=0]=1
            # Only return positive values
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


    def Sr(self, Te=0, Ti=0, Tm=0):
        ''' Temperature-dependent reactant sink 
        try:
            ret = self.Vr - self.Vp - eval(self.K)
        except:
            ret = self.Vr - self.Vp - self.K
        return ret'''
        return self.Vr - self.Vp - eval(self.K)
    
    def decaySr(self, Te=0, Ti=0, Tm=0):
        ''' Temperature-dependent fragmentation 
        # TODO: IS SIGN MISMATCH WARRANTED??!
        try:
            ret = -eval(self.K)
        except:
            ret = self.K
        return ret'''
        return -eval(self.K)

    def depSe(self, Te, *args):
        ''' Electron energy gain '''
        return Te

    def get_nsum(self, ne, ni):
        return max(1,self.e*ne + self.p*ni)
    
    def get_n(self, ne, ni):
#        if (self.e is True) and (self.p is True):
#            print(self.name)
#        return max(self.e*ne * self.p*ni,1)
        if self.e is True: # Electron impact: assume to use for rec. as well
            return ne
        elif self.p is True: # Proton impact
            return ni
        else: # Neither electron or ion impact: assume spontaneous
            return 1

    def interp1d_loglog(self,T):
        from numpy import log10
        return 10**self.interp(log10(T))
        
        
    
    def print_reaction(self):
        ''' Returns formatted a string with the reaction '''
        ret = '{}_{}_{}: '.format(self.database, self.type, self.name) # Append reaction ID
        for i in range(len(self.educts)):    # Loop through all reactants
            ret += (self.e_mult[i] != 1)*'{}*'.format(self.e_mult[i])\
                    + '{}'.format(self.educts[i])\
                    + (i+1 != len(self.educts))*'+ '
        ret += ' => ' # Add reactants to string
        for i in range(len(self.products)):    # Loop through all fragmetns
            ret += (self.p_mult[i] != 1)*'{}*'.format(self.p_mult[i])\
                    + '{}'.format(self.products[i])\
                    + (i+1 != len(self.products))*'+ ' 
            # Add products to string
        
        return ret 


    def pick_rate(self):
        ''' Initialization of the rate attribute '''
        # TODO: Dynamically choose between ni and ne in reaction
        from scipy.interpolate import interp1d
        from numpy import pi, log10
        if self.database == 'ADAS': 
            # TODO: verify ADAS
            print('TBD')
            '''
            if ID == 'EXCITATION': # Electron impact excitation
                # Excitation only possible between current state and nmax
                rn = range(x + 1, self.nmax + 1)   
                fit = 'ADAS' # ADAS-type fit
                Tarr = rdata[database]['T']  # Temperature array for interpolation
            elif ID=='RELAXATION': # Radiative relaxation
                rn = range(1, x)   # Relaxation only possible to lower states
                fit = 'COEFFICIENT' # Coefficient-like decay
                Tarr = None   # Not T-dependent
            # Loop through each of the available final states
            for y in rn:
                try:
                    _name = '{}_{}-{}'.format(ID, x, y)
                    #Get coefficients
                    _ratecoeff = rdata['ADAS']['{}-{}'.format(x, y)] 
                    self.reactions['ADAS'][_name] = Reaction('ADAS',
                            _ratecoeff, fit, data, bg, species, 
                            self.isotope, self.mass, Tarr)
                except:
                    pass
                # Create 1D interpolation function
                self.interpolation = interp1d(
                        self.Tarr, self.coeffs, kind='slinear')   
                return self.ADAS_rate'''
        elif self.database == 'UE': 
            # TODO: verify UE
            print('TBD')
#            t = [i for i in range(self.coeffs.shape[0])]
#            n = [i for i in range(self.coeffs.shape[1])]
#            # Create 2D interpolation function
#            self.interpolation = interp2d(n, t, self.coeffs, kind='linear') 
#            return self.UE_rate
        elif self.database == 'APID': 
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
        elif self.database == 'JOHNSON':
            ''' L.C.JOHNSON, ASTROPHYS. J. 174, 227 (1972). '''
            # TODO: verify against Eirene's COLRAD!
            (i, f) = self.coeffs
            def g(i, f):
                g = [
                        1.133*(f == 1) + 1.0785*(f == 2) + (0.9935 + 0.2328/f
                        - 0.1296/f**2)*(f > 2),
                        -0.4059*(f == 1) - 0.2319*(f==2) - ((0.6282 - 0.5598/f
                        + 0.5299/f**2)/f)*(f > 2),
                        0.07014*(f == 1) + 0.02947*(f == 2) + ((0.3887 - 1.181/f
                        + 1.470/f**2)/f**2)*(f > 2) ]
                x=1 - (f/i)**2
                return g[0] + g[1]/x + g[2]/x**2

            h = 1.054571596e-27
            c = 2.99792458e10
            me = 9.10938188e-28
            e = 4.80274e-10
            I = (me*e**4)/(2*h**2)
            res = (2**6 * e**2 * I**2)/(3**1.5*pi*me*c**3 * h**2)
            freq = (1/f**2 - 1/i**2)
            Afac = (res*g(i, f))/(freq*(i**5)*(f**3))
            self.coeffs = Afac
            return self.coeff_constant

        elif self.database == 'MCCCDB':
            self.interp = interp1d(log10(self.coeffs[:,0]),
                                log10(self.coeffs[:,1]))
            self.xsec = self.interp1d_loglog
            return self.integrated

        elif 'H.' in self.type.upper(): 
            if self.type.upper() == 'H.0':
                print('Potential: TBD')
                return
            elif self.type.upper() == 'H.1':
                self.xsec = self.polyfit_H1
                return self.integrated
            elif self.type.upper() == 'H.2':
                return self.polyfit_H2
            elif self.type.upper() == 'H.3':
                return self.polyfit_H3
            elif self.type.upper() == 'H.4':
                return self.polyfit_H4
            elif self.type.upper() == 'H.5':
                print('Momentum-weighted rates vs. temperature, not in use')
                return
            elif self.type.upper() == 'H.6':
                print('Momentum-weighted rates vs. temperaturei and energy: TBD')
                return
            elif self.type.upper() == 'H.7':
                print('Momentum-weighted rates vs. temperature and density, not in use')
                return
            elif self.type.upper() == 'H.8':
                print('Energy-weighted rates vs. temperature: TBD')
                return
            elif self.type.upper() == 'H.9':
                print('Energy-weighted rates vs. temperature and energy, not in use')
                return
            elif self.type.upper() == 'H.10':
                print('Energy-weighted rates vs. temperature and density: TBD')
                return
            elif self.type.upper() == 'H.11':
                print('Other data: TBD')
                return
            elif self.type.upper() == 'H.12':
                print('Other data: TBD')
                return
            else:
                print('Unknown fit: {}, {}, {}'.format(
                        self.database,self.type,self.name))
        elif self.type.upper() in ['COEFFICIENT', 'RELAXATION']: 
            if ('e' in self.educts) or ('p' in self.educts):
                return self.coeff_rate
            else:
                return self.coeff_constant
        elif self.type.upper() == 'INTERPOLATION': 
            self.Tl, self.Tu = 10**self.coeffs[0,0], 10**self.coeffs[-1,0]
            self.coeffs = interp1d(self.coeffs[:,0], self.coeffs[:,1])
            return self.loglog_interp
        elif self.type.upper() == 'SIGMA': 
            # TODO: verify SIGMA
            print('Sigma rates to be verified!')
            # NOTE: To be verified
            return self.SAWADASIGMA_rate
        else: 
            print('Unknown type "{} {}"'.format(self.database, self.type))


    def ADAS_rate(self, Te, Ti, omegaj=1, **kwargs):
        ''' Function returning the ADAS rate for T '''
        ep = ((self.e + self.p) ==2)
        T = (Te*self.e + Ti*(self.p - ep))/(self.e + self.p - ep)
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


    def extrapolate_polyfit(self, T, fittype, E, ne, ni, frequency = False):
        from numpy import log, exp, log10
        x = log10(T)
        dx=1e-1
        x1 = 0.5
        y1 = log10(fittype(Te=x1,Ti=x1,E=E,ne=ne,ni=ni))
        dy = log10(fittype(Te=x1+dx,Ti=x1+dx,E=E,ne=ne,ni=ni))
        k = -(y1-dy)/(log10(x1+dx)-log10(x1))
        b = y1 - k*log10(x1)
        return self.get_n(ne, ni)**frequency*(10**(k*x+b))


    def polyfit_H1(self, Te, Ti=0, E=0, ne=0, ni=0, extrapolate=True):
        ''' Function returning the EIRENE rate for T '''
        from numpy import log, exp, log10
        ep = ((self.e + self.p) ==2)
        ret = 0
        if (Te < 0.5) and (extrapolate is True):   
            return self.extrapolate_polyfit(Te, self.polyfit_H1, 0, 0, 0, False)
        # Rate coefficient vs temperature
        for i in range(9):
            ret += self.coeffs[i]*(log(Te)**i)
        return self.scale*selg.get_n(ne,ni)**frequency*exp(ret)

    def polyfit_H2(self, Te, Ti, E=0, ne=0, ni=0, frequency=False, extrapolate = True, **kwargs):
        ''' Function returning the EIRENE rate for T '''
        from numpy import log, exp, log10
        ep = ((self.e + self.p) ==2)
        T = (Te*self.e + Ti*(self.p - ep))/(self.e + self.p - ep)
        ret = 0
        if (T < 0.5) and (extrapolate is True):   
            return self.extrapolate_polyfit(T, self.polyfit_H2, E, ne, ni, frequency)
        # Rate coefficient vs temperature
        for i in range(9):
            ret += self.coeffs[i]*(log(T)**i)    
        return self.scale*self.get_n(ne,ni)**frequency*exp(ret)
    
    def polyfit_H3(self, Te, Ti, E=None, ne=0, ni=0, frequency=False, extrapolate=True, **kwargs):
        ''' Function returning the EIRENE rate for T '''
        from numpy import log, exp, log10
        ep = ((self.e + self.p) ==2)
        T = (Te*self.e + Ti*(self.p - ep))/(self.e + self.p - ep)
        ret = 0
        if (T < 0.5) and (extrapolate is True):   
            return self.scale*self.extrapolate_polyfit(T, self.polyfit_H3, E, ne, ni, frequency)
        # Rate coefficient vs temperature and energy
        for i in range(9):
            for j in range(9):
                ret += self.coeffs[i,j]*(log(T)**i)*(log(E)**j)
        return self.scale*self.get_n(ne,ni)**frequency*exp(ret)

    def polyfit_H4(self, Te, Ti, E=0, ne=None, ni=None, frequency=False, extrapolate=True, **kwargs):
        ''' Function returning the EIRENE rate for T '''
        from numpy import log, exp, log10
        ep = ((self.e + self.p) ==2)
        T = (Te*self.e + Ti*(self.p - ep))/(self.e + self.p - ep)
        ret = 0
        if (T < 0.5) and (extrapolate is True):   
            return self.scale*self.extrapolate_polyfit(T, self.polyfit_H4, E, ne, ni, frequency)
        # Rate coefficient vs temperature and density
        for i in range(9):
            for j in range(9):
                ret += self.coeffs[i,j]*(log(T)**i)*(log(ne*1e-8)**j)
        return self.scale*self.get_n(ne,ni)**frequency*exp(ret)



    def loglog_interp(self, Te, Ti=0, E=None, ne=None, ni=None, frequency=False, **kwargs):
        from numpy import log10
        ep = ((self.e + self.p) ==2)
        T = (Te*self.e + Ti*(self.p - ep))/(self.e + self.p - ep)
        if (T >= self.Tl) and (T <= self.Tu):
            fit = self.coeffs(log10(T))
        elif T < self.Tl:
            fit = self.coeffs(log10(self.Tl))
        elif T > self.Tu:
            fit = self.coeffs(log10(self.Tu))
        return self.scale*(self.get_n(ne,ni)**frequency)*10**fit

    def coeff_constant(self, *args, frequency=False, ne=None, ni=None, **kwargs):
        # NOTE: These are always assumed to be in 1/s 
        return self.coeffs*self.scale
            
    def coeff_rate(self, *args, frequency=False, ne=None, ni=None, **kwargs):
        # NOTE: These are always assumed to be in 1/s 
        return self.scale*self.coeffs*self.get_n(ne,ni)**frequency
            

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
        ep = ((self.e + self.p) ==2)
        T = (Te*self.e + Ti*(self.p - ep))/(self.e + self.p - ep)
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
        return self.scale*(4/sqrt(pi))*sqrt((T*ev)/(2*me))*\
                quad(R, 0, inf, args=(T, self.coeffs))[0]
    

    def APID_rate(self, Te, Ti, ne=None, ni=None, frequency=False, **kwargs):
        ''' Function returning the APID rate for T '''
        from numpy import exp, sqrt, log, pi
        ep = ((self.e + self.p) ==2)
        T = (Te*self.e + Ti*(self.p - ep))/(self.e + self.p - ep)
        ''' The APID-4 rates are taken from Janev's 1993 IAEA paper, the 
            analytic solutions from Stotler's svlib routine used in DEGAS2 '''
        I = 13.6/self.coeffs**2
        n = self.get_n(ne, ni)**frequency

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
                return self.scale*(k == 0)*exp(-p)/p + (p == 0)*(1/(k-1))

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
                    return self.scale*ze1
                else:
                    zek = ze1
                    for i in range(2, k+1):
                        zek = (exp(-p) - p*zek)/(i - 1)
                    return self.scale*zek
            
            else:
                kp = int(p + 0.5)
                if k < kp and k <= 10:
                    if kp > 10: 
                        kp = 10
                    zek = en(kp, p)
                    for i in range(kp - 1, k - 1, -1): 
                        zek=(exp(-p) - i*zek)/p
                    return self.scale*zek
            
                else: return self.scale*en(k, p)
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
            return n*6.692e7*sqrt(T)*zi3*(zi1 + zi2)

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

        return n*ret
            

    '''
    def integrand(self, x, T):
        from numpy import exp, log
        return x*self.xsec(x*T)*exp(-x)

    def integrated(self, Te, Ti, ne=0, ni=0,E=None, frequency=False, **kwargs):
        from scipy.integrate import quad
        from numpy import inf, pi
        T = (self.e*Te+self.p*Ti)
        return self.scale*100*self.get_n(ne,ni)**frequency\
                *(4/pi**0.5)*((1.602e-19*T)/\
                (2*9.1093837e-31))**0.5*\
                quad(self.integrand, 0, inf, args=(T))[0]


    '''

    def integrand(self, x, T):
        from numpy import exp, log
        ret = x*self.xsec(x*T)*exp(-x)
        ret[ret==0] = 1e-99
        return log(ret)

    def integrated(self, Te, Ti, ne=0, ni=0, E=None, frequency=False, **kwargs):
        """ Integration in log-log space using trapezoidal rule """
        from numpy import logspace,exp,log, pi, diff
        from scipy.special import logsumexp

        T = (self.e*Te+self.p*Ti)
        x = logspace(-10,4,1000)
        y = self.integrand(x, T)
        deltas = log(diff(x))
        integrand = exp(-log(2.) + logsumexp([logsumexp(y[:-1]+deltas), logsumexp(y[1:]+deltas)]))

        return self.scale*100*self.get_n(ne,ni)**frequency\
                *(4/pi**0.5)*((1.602e-19*T)/\
                (2*9.1093837e-31))**0.5*integrand
