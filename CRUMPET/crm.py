# CRUMPET collisional-radiative model class reaction.py
# Separated from CRUMPET.py by holm10
# Changelog
# 200205 - Separated from CRUMPET.py #holm10
from CRUMPET.tools import Tools


class Crm(Tools):
    ''' Class containing the CRM setup and central routines

    This class stores all information of the CRM setup, as defined in 
    the input file. The class contains all the central routines for
    solving the system of ODEs or creating a Greenland CRM. Inherits
    the methods contained in the Tools class.

    ...

    Attributes
    ----------
    species : dict
        Each entry is a dictionary containing information on the 
        initial density and potential level of the CRM  species. The 
        dictionary keys are used to identify the species
    reactions : dict with databases and Reaction objects included
        A dict of databases, each database consisting of a dict of
        reaction names/handles, which contains the Reaction object: e.g
        reactions[database][reactionname]
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
    path : str
        Path to the parent folder relative to which subsequent paths
        are defined
    isotope : string
        The handle for the isotope being used. Used to identify 
        atomic and molecular species: molecules are taken to be
        {isotope}2{groundstatemol} and atoms to be 
        {isotope}{groundstateatom}
    mass : float
        The isotope mass in AMU, used for calculations in Crm 
    '''


    def __init__(
            self, species, bg, reactionlist ,verbose, NP, path='.', vmax=14,
            nmax=8, rdata=None, isotope='H', mass=1):
        ''' 
        Parameters
        ----------
        species : dict
            a dictionary of species handles, each containing a dict 
            defining the species initial density 'n' and the potential
            of the species 'V'
        background : dict
            a dictionary of the background species, each containing a 
            dict of the species potential 'V'

        reactionlist : dict with databases and reactions
            A dict of databases, each database consisting of a dict of
            reaction names/handles, which contains the reaction data.
            The reactions[database][reactionname] dict requires three
            keys:
            'reactants' : list of strings
                each element is a handle of the reactants. At least one
                of the species must be a background species, and only
                two reactants are currently supported
            'products' : list of strings
                each element is a handle of the products
            'K' : string
                a string of the net kinetic energy transfer for the 
                background reactants to the products. Presently, all
                energy is assumed to end up as ion/atom kinetic energy.
                Supports calculations in the string, e.g. '2*3+4'
        verbose : boolean
            switch whether to print the CRM setup to the prompt: the 
            same info is written to {Crm.path}/logs regardless
        NP : int
            the number of P-species to be considered in the model: the 
            NP first entries in species (e.g. the input file) are taken
            as the P-species. The remainder are assumed Q-species.
        path : string, optional (default: '.')
            a path to which all subsequent paths are relative to. 
            Relative or an absolute path.
        vmax : int, optional (default: 14)
            number of vibrational levels to be considered for 
            vibrationally dependent processes
        nmax : int, optional (default: 8)
            number of electronic levels to be considered for atomic
            electron transitions
        rdata : RateData object, optional (default: None)
            a RateData object that contains the databases and reactions
            defined in the input file reactions
        isotope : string, optional (default: 'H')
            The handle for the isotope being used. Used to identify 
            atomic and molecular species: molecules are taken to be
            {isotope}2{groundstatemol} and atoms to be 
            {isotope}{groundstateatom}
        mass : float, optional (default: 1)
            The isotope mass in AMU, used for calculations in Crm 
        '''
        from os import mkdir,getcwd
        from datetime import datetime

        # Store class objects
        self.species = species
        self.bg = bg
        self.isotope = isotope
        self.mass = mass
        self.slist = list(self.species)
        self.verbose = verbose
        self.Np = NP
        self.path = path
        self.reactions = {} 
        self.vmax = vmax
        self.nmax = nmax
        self.setup_reactions(reactionlist ,rdata)

        # Define reactions for UEDGE raditation
        #self.ionizrad=Reaction('IONIZRAD','UE',
                #self.get_coeff('UE','IONIZRAD'),'UE',['','',''],bg,species,
                #['','',None,None,[0,0,0,0]])
        #self.recrad=Reaction('RECRAD','UE',self.get_coeff('UE','RECRAD'),
                #'UE',['','',''],bg,species,['','',None,None,[0,0,0,0]])
        # Ensure that there is a logs directory under the run path
        try:
            mkdir('{}/logs'.format(self.path))
        except:
            pass
        # Write a log of the CRM setup path/logs
        with open('{}/logs/setup.log'.format(self.path), 'w') as f:
            f.write('CRUMPET run in {} on {}\n'.format(getcwd(),
                    str(datetime.now())[:-7]))
            f.write('Defined species:\n')
            for i in self.slist:
                f.write('    {}\n'.format(i))
            f.write('Defined reactions:\n')
            for dkey, database in self.reactions.items():
                for rkey, reaction in database.items():
                    f.write('{}\n'.format(reaction.print_reaction(dkey, rkey)))
        # Output to stdout if run verbosely
        if self.verbose:
            with open('logs/setup.log', 'rt') as f:
                for m in f:
                    print(m.strip())
        # Do the same for a Diagnostic rate matrix displaying 
        # reaction correlations
        self.DIAGNOSTIC()

# %%%%%%%%%%%%%% INPUT FILE READING TOOLS %%%%%%%%%%%%%%%%%
    def setup_ADAS(self, x, ID, rdata, data):
        ''' Adds an ADAS reactions up to self.nmax'''
        from CRUMPET.reactions import Reaction

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


    def setup_APID(self, x, ID, rdata):
        ''' Adds an APID-style ionization reaction '''
        from CRUMPET.reactions import Reaction
        self.reactions['APID'][self.XY2num(ID,x)] = Reaction('APID', x, 'APID',
                rdata, self.bg, self.species, self.isotope, self.mass)


    def setup_Johnson(self, i, ID, r):
        ''' Adds a Johnson-style relaxation reaction from state n=i'''
        from CRUMPET.reactions import Reaction
        from numpy import pi

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
        for f in range(1, i):
            buff = r.copy()
            res = (2**6 * e**2 * I**2)/(3**1.5*pi*me*c**3 * h**2)
            freq = (1/f**2 - 1/i**2)
            Afac = (res*g(i, f))/(freq*(i**5)*(f**3))
            for label in ['reactants','fragments']:
                buff[label] =  [self.XY2num(a, i, f) for a in buff[label]]
            self.reactions['JOHNSON'][self.XY2num(ID, i, f)] = Reaction(
                    'JOHNSON', Afac, 'COEFFICIENT', buff, self.bg, 
                    self.species, self.isotope, self.mass)


    def setup_reactions(self, reactionlist, rdata):
        ''' Adds the reactions to the CRM

        Called by __init__, adds the reactions from reactionlist to
        self.reactions using the rates in rdata

        Parameters
        ----------
        reactionlist : list 
            a list of lists for setting up the CRM reactions. Each 
            element is a list of the following data
            handle : string
                database_reaction as specified in the input
            reactants : list of strings
                each element is a handle of the reactants. At least one
                of the species must be a background species, and only
                two reactants are currently supported
            products : list of strings
                each element is a handle of the products
            K : string
                a string of the net kinetic energy transfer for the 
                background reactants to the products. Presently, all
                energy is assumed to end up as ion/atom kinetic energy.
                Supports calculations in the string, e.g. '2*3+4'
        rdata : RateData object, optional (default: None)
            a RateData object that contains the databases and reactions
            defined in the input file reactions

        Returns
        -------
        None
        '''
        from CRUMPET.reactions import Reaction
        from numpy import zeros, pi
        
        ''' LOOP OVER ALL THE DEFINED ReactionS '''
        for database, reaction in reactionlist.items():
            # Create database if not already created
            try:
                self.reactions[database]
            except:
                self.reactions[database] = {}
            for ID, data in reaction.items():
                # Split the definitions into names and databases
                ndep = ('N=$' in ''.join(data['reactants'] + 
                        data['fragments']).upper())
                vdep = ('V=$' in ''.join(data['reactants'] + 
                        data['fragments']).upper())
                # Loop through states, if necessary. Dynamicallt set boundaries
                # according to electronic or vibrational transitions
                for x in range(ndep, 1 + vdep*self.vmax + ndep*self.nmax):
                    # Vibrational/electronic dependence present
                    if '$' in ID:
                        buff = data.copy()
                        # Substitute state into reactants and product strings 
                        for label in ['reactants', 'fragments']:
                            buff[label] = [self.XY2num(i, x) for i in 
                                    buff[label]]
                        # %%% ADAS rates detected %%%
                        if database == 'ADAS':
                            self.setup_ADAS(x, ID, rdata, rlist, plist, 
                                    data['K'])
                        # %%% APID rates detected %%%
                        elif database == 'APID':
                            self.setup_APID(x, ID, buff)
                        # %%% APID rates detected %%%
                        elif database == 'JOHNSON':
                            self.setup_Johnson(x, ID, buff)
                        # %%% Neither of the above: rate for transitions %%%
                        else:
                            # Assume ladder-like vibrational transitions (+/-1)
                            for y in range(-1*('&' in ID), 2, 2):
                                # Limit transitions to [0,vmax]
                                if x + y in range(self.vmax + 1):
                                    rebuff = buff.copy()
                                    # Retain intial and final states in name
                                    vID = self.XY2num(ID, x, x + y) 
                                    for label in ['reactants', 'fragments']:
                                        rebuff[label] = [self.XY2num(i, x, 
                                                x + y) for i in rebuff[label]]
                                    _name = vID
                                    _database = database
                                    _ratecoeff = rdata[database][vID] 
                                    _rtype = 'RATE'
                                    _rlist = rebuff
                                    self.reactions[_database][_name] = \
                                            Reaction(_database, _ratecoeff, 
                                            _rtype, _rlist, self.bg, 
                                            self.species, self.isotope, 
                                            self.mass)
                    #%%%%% Read custom rates %%%%%
                    elif database == 'CUSTOM':
                        self.setup_custom(ID.strip(), database)
                    #%%% EIRENE/UEDGE-DEGAS rates %%%
                    elif database in ['HYDHEL', 'AMJUEL', 'H2VIBR', 'UE']:
                        self.reactions[database][ID] = Reaction(database,
                                rdata[database][ID], 'RATE'*(database != 'UE') + 
                                'UE'*(database == 'UE'), data, self.bg, 
                                self.species, self.isotope, self.mass)
                    #%%% Fell through loop %%%
                    else:
                        print(  'Database "{}" not recognized!\n'
                                'Aborting.'.format(database))
                        return


    def setup_custom(self, fname, database):
        ''' Adds the custom reactions defined in the input to the CRM
        
        Parses any custom input rate file defined in the input file
        and adds their Reaction objects to the CRM.
    
        Parameters
        ----------
        fname : string
            path to the custom rate setup file relative to Crm.path
        database : string
            database handle for the rates in the custom rate file
        
        Returns
        -------
        None
        '''
        from CRUMPET.reactions import Reaction
        from numpy import zeros, pi
        # Parse the custom rate file into a list and retain subcards
        data, _, subcards = self.file2list(self.path, fname)
        subcards = subcards[1:]
        _database = data[0].split()[1].strip() # Database is defined as card
        try:
            self.reactions[_database]
        except:
            self.reactions[_database] = {}
        # Loop through the rates, each a separate subcard
        for i in range(len(subcards) - 1):
            subc = subcards[i] # Helper index
            # Rate typy, first subcard entry
            fit = data[subc].split()[1].upper() 
            name = data[subc].split()[2] # Name/ID of reacrtion - second entry
            rdata = {
                    'reactants': data[subc + 1].split(' > ')[0],
                    'fragments': data[subc + 1].split(' > ')[1],
                    'K': '0',
                    }
            # Execute if no vib dependence, loop if vibr. dep. process
            for j in range(1 + ('$' in name)*self.vmax):
                buff = rdata.copy()
                # %%% Vibrationally dependent process %%%
                if '$' in name:
                    # Read kinetic energy for each process
                    for m in range(0, 100):
                        if data[m][0].upper() == 'K': 
                            buff['K'] = data[m].strip().split('=')[-1]
                        elif data[m][0].upper() == 'V':
                            break
                    # Write data
                    _name = self.XY2num(name, j)
                    for label in ['reactants', 'fragments']:
                        buff[label] = self.XY2num(buff[label], 
                                j).strip().split(' + ') 
                    # Coefficients for v-level
                    _ratecoeff = [float(x) for x in 
                            data[subc + m+  j*2].split()] 
                # %%% Specified rate %%%
                else:
                    # Read the kinetic energy of the process
                    for m in range(2, 100):
                        if data[subc + m][0] == 'K': 
                            buff['K'] = data[subc + m].strip().split('=')[-1]
                        else:
                            break
                    # Cross-section as defined in SAWADA 95 has special form
                    if fit == 'SIGMA': 
                        _ratecoeff = [float(x) for x in data[subc + m].split()] 
                    # Other processes have pre-defined form
                    else: 
                        _ratecoeff = float(data[subc + m])
                    # Write data
                    _name = name
                    for label in ['reactants', 'fragments']:
                        buff[label] = buff[label].strip().split(' + ') 
                # Store reaction
                self.reactions[_database][_name] = Reaction(_database,
                        _ratecoeff, fit, buff, self.bg, self.species,
                        self.isotope, self.mass)

# %%%%%%%%%%%%%% ENF OF INPUT FILE READING TOOLS %%%%%%%%%%%%%%%%%

    def populate(self, mode, Te, ne, Ti=None, ni=None, E=0.1, Tm=0):
        ''' Function populating and returning the rate matrix and source

        Parameters
        ----------
        mode : string
            Defines the type of matrix to be written according to
            'diagnostic' - 2D list of reaction handles
            'R' - rate coefficient matrix in cm**3 s**-1
            'M' - rate matrix in s**-1
            'Sgl' - energy rate matrix in eV s**-1
            'I' -   gamma count matrix, marking gammas s**-1 for
                    off-diagonal transitions  
            'E' -   gamma energy matrix, marking the gamma energies
                    for the off-diagonal transitions in eV
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        Ti : float, optional (default: None)
            ion temperature in eV. Default sets Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Default sets ni=ne
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        Tm : float, optional (default: 0)
            local molecular temperature for calculating energy losses 
            associated with molecules. Defaults to E if Tm=0.

        Returns
        -------
        mat : ndarray
            Rate matrix containing the rates at the given 
            temperature and density. A 2D list in case mode=diagnostic
        ext : ndarray
            Vector of external sinks/sources. List if mode=diagnostic
        '''
        from numpy import zeros,array,sum,transpose
        
        N = len(self.species)
        if mode == 'diagnostic':
            # Setup a 2D diagnostic list for the matrix and a list for 
            # the external source
            ext_source = [[] for x in range(N)]
            ret = [[[] for y in range(N)] for x in range(N)]
        elif mode in ['E', 'I', 'Sgl']:
            ret = zeros((N, N, 2 + 3*(mode == 'Sgl')))
            ext_source = zeros((N, 2+3*(mode == 'Sgl')))
        else:
            # Setup a matrix and vector
            ret = zeros((N, N))
            ext_source = zeros(N)
        for i in range(len(self.species)):
            #%%% Walk through each row (species) %%%
            for dkey, database in self.reactions.items():
                for rkey, r in database.items():
                    #%%% Sort the species of each reaction into %%%
                    #%%%  the appropriate column %%%
                    Sgl = r.getS(Te, Ti, Tm, E)
                    # TODO: what if three-particle reaction?
                    bg = r.e*ne + r.p*ni # Specify density for reactions
                    ext = 0
                    if mode != 'diagnostic':
                        bgm = (r.e*ne) * (r.p*ni)
                    if mode == 'R':
                        bgm = r.rate(Te, Ti, E=E, ne=ne)
                        bg = r.rate(Te, Ti, E=E, ne=ne)
                    elif mode in ['M']:#,'Sgl']:
                        # Specify density for external source
                        bgm = bgm*r.rate(Te, Ti, E=E, ne=ne) 
                        # Assure that auto-processes are considered
                        bg = max(bg, 1)*r.rate(Te, Ti, E=E, ne=ne)    
                    elif mode == 'Sgl':
                        # Specify density for external source
                        bgm = bgm*r.rate(Te, Ti, E=E, ne=ne)*Sgl[:,0] 
                        # Assure that auto-processes are considered
                        bg = max(bg, 1)*r.rate(Te, Ti, E=E, ne=ne)*Sgl[:,0]
                        ext = Sgl[:,1]
                    elif mode == 'I':
                        # Specify density for external source
                        bgm = bgm*r.rate(Te, Ti, E=E, ne=ne)*Sgl[-2:,0] 
                        # Assure that auto-processes are considered
                        bg = max(bg, 1)*r.rate(Te, Ti, E=E, ne=ne)*\
                                (abs(Sgl[-2:,0]) > 0)    
                    elif mode == 'E':
                        bgm = Sgl[-2:,0] # Specify density for external source
                        # Assure that auto-processes are considered
                        bg = Sgl[-2:,0]  
                    j = None # Set flag to identify external sources
                    # Loop through each reaction defined in the CRM
                    for rea in range(len(r.reactants)):
                        # Find the column into which the fragments goes: 
                        # if background mark external
                        try:
                            # Get the product species index of the 
                            # correct column
                            j = self.slist.index(r.reactants[rea])  
                        except:
                            continue
                        # If the species (row index) is a reactant, r is a 
                        # depletion process 
                        if self.slist[i] == r.reactants[rea]:   
                            ''' DEPLETION '''
                            if mode == 'diagnostic':
                                ''' Diagnostic matrix '''
                                # Print the rate to the correct element
                                ret[i][i].append('-{}*{}_{}{}'.format(
                                        r.r_mult[rea], dkey, rkey, bg))
                            elif mode in ['R','M']:
                                ''' Rate coefficient matrix '''
                                # Calculate the Rate coefficient and store 
                                # appropriately
                                ret[i,i]-=r.r_mult[rea]*bg 
                    # Loop through the reaction fragments
                    for frag in range(len(r.fragments)):    
                        ''' SOURCE '''
                        # Do nothing if background fragment
                     #   if r.fragments[frag] not in self.slist: 
                     #       continue
                        # If fragment enters row, add in appropriate column as
                        # defined by reactant
                        if (r.fragments[frag] in self.slist) and (
                                self.slist.index(r.fragments[frag]) == i):
                            multiplier = r.f_mult[frag]**(mode not in 
                                    ['Sgl', 'I', 'E'])   # Fragment multiplier
                            # External flag triggered, store to external source
                            if j is None: 
                                ''' EXTERNAL SOURCE '''
                                if mode == 'diagnostic':
                                    ''' Diagnostic matrix '''
                                    ext_source[i].append('+{}*{}_{}{}'.format(
                                            multiplier, dkey, rkey, bg))
                                else:
                                    try:
                                        ext_source[i] += multiplier*bgm + ext
                                    except:
                                        ext_source[i,:] += multiplier*bgm + ext
                            # No trigger of external source, store to 
                            # appropriate location in matrix
                            else: 
                                ''' INTERNAL SOURCE '''
                                if mode == 'diagnostic':
                                    ''' Diagnostic matrix '''
                                    ret[i][i].append('+{}*{}_{}{}'.format(
                                            multiplier, dkey, rkey, bg))
                                else:
                                    try:
                                        ret[i,j] += multiplier*bg + ext
                                    except:
                                        ret[i,j,:] += multiplier*bg + ext
        return ret, ext_source

        
    def write_matrix(self, mat, ext, char, te, ne, ti, ni, E, form='{:1.1E}'):
        ''' Writes the matrix to the logs
        
        Parameters
        ----------
        mat : ndarray or 2D list
            square matrix or list with data to be written
        ext : ndarray or list
            1D vector with external data to be written
        char : string
            Identifier for file to be written
        te : float
            Electron temperature for header
        ne : float
            Electron density for header
        ti : float
            ion temperature for header
        ni : float
            ion density for header
        E : float
            molecule energy for header
        form : string, optional (default: ':1.1E')
            formatting for numerical output

        Returns
        -------
        None
        '''
        from os import getcwd
        from datetime import datetime

        Slabels = ['S_e', 'S_ia', 'S_V', 'S_g']
        # Create a file to write to according to char
        with open('{}/logs/{}_matrix.log'.format(self.path,char), 'w') as f:  
            if char == 'R': # Rate coefficient  matrix is being written
                f.write('Diagnostic rate coefficient (density-independend) '
                        'matrix for CRUMPET run in {} on {}\n'.format(getcwd(),
                        str(datetime.now())[:-7]))
                f.write('Te={} eV, Ti={} eV, E={} eV\n'.format(te,ti,E))
            elif char == 'M': # Rate matrix is being written
                f.write('Diagnostic rate (density-dependend) matrix for CRUMPET '
                        'run in {} on {}\n'.format(getcwd(),
                        str(datetime.now())[:-7]))
                f.write('Te={} eV, Ti={} eV, ne={} 1/cm**3, ni={} '
                        '1/cm**3, E={} eV\n'.format(te,ti,ne,ni,E))
            elif char == 'S': # Rate matrix is being written
                f.write('Diagnostic energy loss (density-dependend) matrix '
                        'for CRUMPET run in {} on {}\n'.format(getcwd(),
                        str(datetime.now())[:-7]))
                f.write('Te={} eV, Ti={} eV, ne={} 1/cm**3, ni={} 1/cm**3, '
                        'E={} eV\n'.format(te, ti, ne, ni, E))
            # Create header line
            out = '{}-MAT|'.format(char.upper()).rjust(10) 
            for s in self.slist:
                out += s.rjust(10,' ')
            f.write(out + '\n'+'_'*(1 + len(self.species))*10 + '\n')

            # Loop through the matrix
            for m in range(len(mat)):
                if len(mat) == len(self.species):
                    # Create a tabulated file with the row species displayed
                    out = self.slist[m]+'|' 
                else:
                    # Create a tabulated file with the row species displayed
                    out = Slabels[m]+'|' 
                out = out.rjust(10,' ')
                # Add each element to the line with 10 characters reserved
                for e in mat[m,:]:
                    if e == 0:
                        out += ' '*9+'_'
                    else:
                        out += form.format(e).rjust(10,' ')
                f.write(out + '\n')
            # Write the external source to the bottom of the output
            f.write('_'*(1 + len(self.species))*10 + '\n')
            out = 'S_ext|'.rjust(10,' ')
            for s in ext:
                out += form.format(s).rjust(10,' ')
            f.write(out)


    def DIAGNOSTIC(self):
        ''' Writes a log of the diagnostic reaction matrix '''
        from os import getcwd
        from datetime import datetime
 
        # Get the 2D diagnostic list and external source list
        dia, ext = self.populate('diagnostic', 0, '*[ne]', 0, '*[ni]') 
        # Write the diagnostic matrix to the logs
        with open('{}/logs/reaction_matrix.log'.format(self.path), 'w') as f:
            # Write the header line
            f.write('Diagnostic reaction matrix for CRUMPET run in {} on '
                    '{}\n'.format(getcwd(), str(datetime.now())[:-7]))
            # Loop through each species
            for i in range(len(dia)):
                # Write the species header
                f.write('\n\n\n======== {} ========\n'.format(self.slist[i]))
                # Start with the sinks
                f.write('{} DEPLETION:\n                 {}\n'.format(
                        self.slist[i],dia[i][i]))
                # Then do the sources
                for j in range(len(dia[i])):
                    if len(dia[i][j])>0:    # Don't write anything empty
                        if i != j:    # Don't duplicate depletion
                            f.write('From {}:\n                 {}\n'.format(
                                    self.slist[j], dia[i][j]))
                # Write external source last, if applicable
                if len(ext[i]) > 0:   
                    f.write('{} EXTERNAL SOURCE:\n'
                            '                 {}\n'.format(self.slist[i],
                            ext[i]))
        # If verbose, output the log content
        if self.verbose:
            with open('logs/reaction_matrix.log', 'rt') as f:
                for m in f:
                    print(m.strip())
            

    def R(self, Te, Ti=None, E=0.1, sparse=False, write=True):
        ''' Returns the rate coefficients in cm**3 s**-1

        Parameters
        ----------
        Te : float
            electron temperature in eV
        Ti : float, optional (default: None)
            ion temperature in eV. Defaults to Ti=Te
        E : float, optional (default: 0.1)
            molecular energy in eV for rates
        sparse : boolean, optional (default: False)
            Switch to use sparse arrays (True) or not (False)
        write : boolean, optional (default: True)
            Switch for writing the R-matrix to the logs

        Returns
        -------
        mat, ext
        mat : ndarray
            Rate matrix containing the rates at the given temperature
            and density in cm**3 s**-1
        ext : ndarray
            Vector of external sinks/sources in cm**3 s**-1
        '''
        from scipy.sparse import csc_matrix

        if Ti is None: 
            Ti=Te # Check for Ti, set if necessary
        # Set densities to one to get rate coefficients as output
        ni=1
        ne=1
        R, ext = self.populate('R', Te, 0, Ti, 0, E)
        if write: # Write to log if requested
            self.write_matrix(R, ext, 'R', Te, 0, Ti, 0, E)
            if self.verbose: # Print rate matrix to stdout if running verbose
                with open('logs/R_matrix.log', 'rt') as f:
                    for m in f:
                        print(m.strip())
        if sparse: 
            R = csc_matrix(R)  # Use sparse format if requested
        return R, ext


    def M(self, Te, ne, Ti=None, ni=None, E=0.1, sparse=False, write=True):
        ''' Returns the rate coefficients in cm**3 s**-1

        Parameters
        ----------
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**3
        Ti : float, optional (default: None)
            ion temperature in eV. Defaults to Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Defaults to ni=ne
        E : float, optional (default: 0.1)
            molecular energy in eV for rates
        sparse : boolean, optional (default: False)
            Switch to use sparse arrays (True) or not (False)
        write : boolean, optional (default: True)
            Switch for writing the R-matrix to the logs
            
        Returns
        -------
        mat, ext
        mat : ndarray
            Rate matrix containing the rates at the given temperature
            and density in cm**3 s**-1
        ext : ndarray
            Vector of external sinks/sources in s**-1
        '''
        from scipy.sparse import csc_matrix

        if Ti is None: 
            Ti = Te # Check for Ti, set if necessary
        if ni is None: 
            ni = ne # Check for ni, set if necessary
        M, ext = self.populate('M', Te, ne, Ti, ni, E)
        if write:   # Write to log if requested
            self.write_matrix(M, ext, 'M', Te, ne, Ti, ne, E)
            if self.verbose: # Print output if running verbose
                with open('logs/M_matrix.log', 'rt') as f:
                    for m in f:
                        print(m.strip())
        if sparse: 
            M = csc_matrix(M) # Use sparse format if requested
        return M,ext
 

    def Sgl(self, Te, ne, Ti=None, ni=None, E=0.1, Tm=0, write=False):
        ''' Returns the energy rates in eV s**-1

        Parameters
        ----------
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**3
        Ti : float, optional (default: None)
            ion temperature in eV. Defaults to Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Defaults to ni=ne
        E : float, optional (default: 0.1)
            molecular energy in eV for rates
        sparse : boolean, optional (default: False)
            Switch to use sparse arrays (True) or not (False)
        write : boolean, optional (default: True)
            Switch for writing the R-matrix to the logs
            
        Returns
        -------
        mat, ext
        mat : ndarray
            Rate matrix containing the energy rates at the given 
            temperature and density in eV s**-1
        ext : ndarray
            Vector of external energy sinks/sources in eV
        '''
        from numpy import matmul, block
        from numpy.linalg import inv

        if Ti is None: 
            Ti=Te # Check for Ti, set if necessary
        if ni is None: 
            ni=ne # Check for ni, set if necessary
        mat, ext = self.populate('Sgl', Te, ne, Ti, ni, E, Tm=Tm)    
        if write:
            title = ['Sgl_el', 'Sgl_ia', 'Sgl_v', 'Sgl_ga', 'Sgl_gm']
            for i in range(5):
                self.write_matrix(mat[:,:,i], ext[:,i], title[i], Te, ne, Ti, 
                        ni, E, form='{:1.2E}')
        U = [
                [mat[:,:,0], ext[:,0]],
                [mat[:,:,1], ext[:,1]],
                [mat[:,:,2], ext[:,2]],
                [mat[:,:,3], ext[:,3]],
                [mat[:,:,4], ext[:,4]], 
            ]
        # Get the full rate matrix
        M, G = self.M(Te, ne, Ti, ni, E, write=write) 
        MP = M[:self.Np, :self.Np]
        MQ = M[self.Np:, self.Np:]
        V = M[self.Np:, :self.Np]
        H = M[:self.Np, self.Np:]
        ret = []
        for S in U:
            UP = S[0][:self.Np, :self.Np]
            UQ = S[0][self.Np:, self.Np:]
            UV = S[0][self.Np:, :self.Np]
            UH = S[0][:self.Np, self.Np:]
            P = UP - matmul(UH, matmul(inv(MQ), V))
            Q = UV - matmul(UQ, matmul(inv(MQ), V))
            ret.append([block([[P], [Q]]), S[1]])
        return ret


    def I(self, Te, ne, Ti=None, ni=None, E=0.1, Tm=0, write=False):
        ''' Returns the photon intensity in s**-1

        Parameters
        ----------
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**3
        Ti : float, optional (default: None)
            ion temperature in eV. Defaults to Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Defaults to ni=ne
        E : float, optional (default: 0.1)
            molecular energy in eV for rates
        sparse : boolean, optional (default: False)
            Switch to use sparse arrays (True) or not (False)
        write : boolean, optional (default: True)
            Switch for writing the R-matrix to the logs
            
        Returns
        -------
        # TODO: remove external sources? No menaing
        [ [mata, exta], [matm, extm]]
        mata : ndarray
            Intensity of atomic lines between the row- and column
            species
        exta : ndarray
            Intensity of atomic lines from external source
        matm : ndarray
            Intensity of molecular lines between the row- and column
            species
        extm : ndarray
            Intensity of molecular lines from external source
        '''
        if Ti is None: 
            Ti=Te # Check for Ti, set if necessary
        if ni is None: 
            ni=ne # Check for ni, set if necessary
        mat, ext = self.populate('I', Te, ne, Ti, ni, E, Tm=Tm)    
        if write:
            title=['Ia', 'Im']
            for i in range(5):
                self.write_matrix(mat[:,:,i], ext[:,i], title[i], Te, ne, Ti,
                        ni, E, form='{:1.2E}')
        return  [[mat[:,:,0], ext[:,0]], [mat[:,:,1], ext[:,1]]]


    def E(self, Te, ne, Ti=None, ni=None, E=0.1, Tm=0, write=False):
        ''' Returns the photon power in eV s**-1

        Parameters
        ----------
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**3
        Ti : float, optional (default: None)
            ion temperature in eV. Defaults to Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Defaults to ni=ne
        E : float, optional (default: 0.1)
            molecular energy in eV for rates
        sparse : boolean, optional (default: False)
            Switch to use sparse arrays (True) or not (False)
        write : boolean, optional (default: True)
            Switch for writing the R-matrix to the logs
            
        Returns
        -------
        # TODO: remove external sources? No menaing
        [ [mata, exta], [matm, extm]]
        mata : ndarray
            Power of atomic lines between the row- and column
            species
        exta : ndarray
            Power of atomic lines from external source
        matm : ndarray
            Power of molecular lines between the row- and column
            species
        extm : ndarray
            Power of molecular lines from external source
        '''
        if Ti is None: 
            Ti=Te # Check for Ti, set if necessary
        if ni is None: 
            ni=ne # Check for ni, set if necessary
        mat, ext = self.populate('E', Te, ne, Ti, ni, E, Tm=Tm)    
        if write:
            title=['Ea', 'Em']
            for i in range(5):
                self.write_matrix(mat[:,:,i], ext[:,i], title[i], Te, ne, Ti,
                        ni, E, form='{:1.2E}')
        return  [[mat[:,:,0], ext[:,0]], [mat[:,:,1], ext[:,1]]]

 
    def solve_crm(self, t, Te, ne, Ti=None, ni=None, E=0.1, Tm=0, Sext=True, gl=True,
            n=None, Qres=True, densonly=False, write=False):
        ''' Solves the time-evolution of the CRM system at given T and n

        Parameters
        ----------
        t : float
            end-time of simulation in s
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        Ti : float, optional (default: None)
            ion temperature in eV. Default sets Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Default sets ni=ne
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        Sext : boolean, optional (default: True)
            switch for considering 'external' sinks/sources (e.g. 
            re-association)
        Tm : float, optional (default: 0)
            local molecular temperature for calculating energy losses 
            associated with molecules. Defaults to E if Tm=0.
        gl : boolean, optional (default: False)
            switch determining wheter to evaluate the Np x Np 
            Greenland CRM (True) or solve the full NxN system of 
            ODEs (False)
        n0 : list of floats, optional (default: None)
            Initial density distribution of the species defined in the 
            input. List must be of the same length as the species list,
            and the order is determined by the input file. Defaults to 
            the initial densities defined in the input if n0 is None.
        Qres : boolean, optional (default: True)
            Switch determining whether to resolve the Q-species energy
            transfer or return the net energy transfer divided into 
            processes. 
        densonly : boolean, optional (default: False)
            Only returns density rates if True, otherwise evaluates
            both density and energy rates
        
        Returns
        -------
        scipy.integrate._ivp.ivp.OdeResult object
            The solution of the IVP ODE problem up to time t. A function 
            that interpolates the solution in the interval [0,t] is
            accessed by the method sol().
            The solution dimensionality along axis=0 is determined by 
            the CRM setup and along axis=1 by the time-resolution 
            according to the following:
                density=True 
                    Qres=False : Np
                        Data ordered according to input species list
                    Qres=True : N (requires gl=False)
                        Data ordered according to input species list

                density=False
                    Qres=False
                        gl=True : 5 + Np
                            Total Ee, Ei/a, Epot, Erada, Eradm, followed
                            by Np-resolved density rates
                        gl=False : 5 + N
                            Total Ee, Ei/a, Epot, Erada, Eradm, followed
                            by species-resolved density rates
                    Qres=True
                        gl=True : 5*N + Np
                            Species-resolved Ee, Ei/a, Epot, Erada, 
                            Eradm, followed by Np-resolved density rates
                        gl=False : 6*N
                            Species-resolved Ee, Ei/a, Epot, Erada, 
                            Eradm, followed by species-resolved density
                            rates
        '''
        from numpy import block, zeros, matmul, reshape, sum
        from  numpy.linalg import inv
        from scipy.integrate import solve_ivp

        N = len(self.species)
        Np = self.Np
        if n is None:
            n0 = self.n0()
        else:
            n0 = n
        # Get the full rate matrix
        M, G = self.M(Te, ne, Ti, ni, E, write=write) 
        if gl is True:
            M, G, n0 = self.gl_crm(M, G, Sext, n0)
            sfunc = self.Sgl
        else:
            sfunc = self.S
        if densonly is False:
            U=sfunc(Te, ne, Ti, ni, E, Tm, write=write)
            if Qres is True: 
                mat = block([
                        [zeros((N, 5*N)), U[0][0]],
                        [zeros((N, 5*N)), U[1][0]],
                        [zeros((N, 5*N)), U[2][0]],
                        [zeros((N, 5*N)), U[3][0]],
                        [zeros((N, 5*N)), U[4][0]],
                        [zeros((Np*gl + N*(not gl), 5*N)), M],
                        ])
                n = zeros((len(mat),))
                n[-Np*gl - N*(not gl):] = n0
                ext = block([U[0][1], U[1][1], U[2][1], U[3][1], U[4][1], G])
            else:
                mat = block([
                        [zeros((5,)), sum(U[0][0], axis=0)],
                        [zeros((5,)), sum(U[1][0], axis=0)],
                        [zeros((5,)), sum(U[2][0], axis=0)],
                        [zeros((5,)), sum(U[3][0], axis=0)],
                        [zeros((5,)), sum(U[4][0], axis=0)],
                        [zeros((Np*gl + N*(not gl), 5)), M],
                        ])
                if gl is True:
                    n = zeros((len(mat),))
                    n[-Np:] = n0
                else:
                    n = block([zeros((5,)), n0])
                    n = n.reshape((N+5,))
                ext = block([
                        sum(U[0][1], axis=0), sum(U[1][1], axis=0), 
                        sum(U[2][1], axis=0), sum(U[3][1], axis=0), 
                        sum(U[4][1],axis=0), G
                        ])
        else:
            mat = M
            ext = G
            n = n0
        return solve_ivp(lambda x, y: self.ddt(y, x, mat, ext), (0, t),
                n, 'LSODA', dense_output=True)

    def old_steady_state(
            self, Te, ne, na, nm, diva, divm, srca, srcm, bga, bgm, psor,
            psorgc, vol, E=0.1, Ti=None, ni=None, Sext=True, gl=False,
            plot=False):
        ''' Solves the steady-state population of atoms and molecules
        Returns the ground state atom, molecule, and total vibrationally
        excited molecule populations. Assumes the reactions defined for 
        in Crm and the external sinks and sources defined and simulates 
        the population for 10s, assuming steady-staten to be achieved by
        that time.
        Parameters
        ----------
        Te : float
            electron temperature in J
        ne : float
            electron density in m**-3
        na : float
            initial guess for atom density in m**-3
        nm : float
            initial guess for molecule density in m**-3
        diva : float
            atom divergence of domain (transport in/out) in s**-1
        divm : float
            molecule divergence of domain (transport in/out) in s**-1
        srca : float
            atom particle source in domain in s**-1
        srcm : float
            molecule particle source in domain in s**-1
        bga : float
            atom background source in domain in s**-1
        bgm : float
            molecule background source in domain in s**-1
        psor : float
            ionization sink in domain in s**-1
        psorgc : float
            recombination source in domain in s**-1
        vol : float
            volume of the domain in m**3
        Ti : float, optional (default: None)
            ion temperature in eV. Default sets Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Default sets ni=ne
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        Sext : boolean, optional (default: True)
            switch for considering 'external' sinks/sources (e.g. 
            re-association)
        gl : boolean, optional (default: False)
            switch determining wheter to evaluate the Np x Np 
            Greenland CRM (True) or solve the full NxN system of ODEs (False)
        plot : boolean, optional (default: False)
            show a plot of the 10s time-evolution
        Returns
        -------
        na, nm, nmvib
        na : float
            steady-state atom density in m**-3
        nm : float
            steady-state atom density in m**-3
        nmvib - float
            total density of all vibrationally excited molecular species 
            in m**-3
        '''
        # TODO: Solve steady-state from dn/dt=0!
        from numpy import log, ones, matmul, polyfit, zeros, exp, sum
        from scipy.optimize import fsolve, minimize
        from scipy.linalg import norm
        from scipy.integrate import solve_ivp
        from matplotlib.pyplot import figure

        def dt(t, n, mat, ext, n0):
            return matmul(mat, n) + ext

        Te /= 1.602e-19
        ne /= 1e6
        # Set volume to be cm**3
        vol *= 1e6
        # Set initial density in [1/cm**3]
        n = zeros((len(self.n0()),))
        n[0] = 1e-6*na
        n[1] = 1e-6*nm
        # Get the rate matrices 
        mat, ext = self.M(Te, ne, Ti, ni, E, write=True) 
        if gl is True:
            # Get matrices and eigenvalues
            mat, ext, n = self.gl_crm(mat, ext, n=n) 
        # Set external sources to be [1/cm**3/s]
        ext[0] = (psorgc + diva + srca + bga)/vol
        ext[1] = (divm + srcm + bgm)/vol
        # Normalize ionzation rate to UEDGE values 
        mat[0,0] = (psor/ne)/vol
        # Simulate to SS
        ss_sol = solve_ivp(lambda x,y: dt(x, y, mat, ext, n), (0, 10), n,
                'LSODA', dense_output=True)
        if plot is True:
            f = figure(figsize = (7,7))
            ax = f.add_subplot(111)
            for i in range(len(n)):
                line = '-' + '-'*('H2' in self.slist[i])
                ax.semilogx(ss_sol.t, ss_sol.y[i,:], line, label=self.slist[i])
                if i == 0:
                    ax.semilogx(ss_sol.t, ss_sol.y[i,:], 'k.', 
                            label=self.slist[i]) 
            ax.legend(ncol=3)
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Density [cm**-3]')
        return [ss_sol.y[0,-1]*1e6, ss_sol.y[1,-1]*1e6, 
                1e6*self.totmol(ss_sol.y[:,-1])]

    def get_relv(self, Te, ne, E=0.1, Ti=None, ni=None):
        sol = self.steady_state(Te, ne, E, Ti, ni)
        return sol/sol[0]

    def steady_state(
            self, Te, ne, E=0.1, Ti=None, ni=None, Sext=True, 
            gl=False, plot=False, n0=None, ext=None, ioniz=0, dt=False):
        ''' TODO: update documentation

        Solves the steady-state population of atoms and molecules

        Returns the ground state atom, molecule, and total vibrationally
        excited molecule populations. Assumes the reactions defined for 
        in Crm and the external sinks and sources specified and solves
        the steady-state system.

        Parameters
        ----------
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        vol : float
            volume of the domain in cm**3
        Ti : float, optional (default: None)
            ion temperature in eV. Default sets Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Default sets ni=ne
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        Sext : boolean, optional (default: True)
            switch for considering 'external' sinks/sources (e.g. 
            re-association)
        gl : boolean, optional (default: False)
            switch determining wheter to evaluate the Np x Np 
            Greenland CRM (True) or solve the full NxN system of ODEs 
            (False)
        plot : boolean, optional (default: False)
            show a plot of the 10s time-evolution
        n0 : ndarray, optional (default: None)
            Initial density distribution of the species defined in the 
            input. List must be of the same length as the species list,
            and the order is determined by the input file. Defaults to 
            the initial densities defined in the input if n0 is None.
        ext : ndarray, optional (default: None)
            External sinks and sources for the species, not considered 
            in CRUMPET, including transport, artificial sinks and 
            sources, etc. Must be of length Nx1, as defined by the input
            species. Defaults to zeros, given in cm**-3 s**-1
        ioniz : float, optional (default: None)
            Additional ionization source (<0) to applied to the first 
            species. Necessary for evaluating CRM-transport code 
            couplings where ionization is handled by the auxilliary
            code, and not included in the CRUMPET model. Given in 
            cm**-3 s**-1

        Returns
        -------
        n_ss : ndarray
            steady-state populations
        '''
        from numpy import log, ones, matmul, polyfit, zeros, exp, sum
        from scipy.optimize import fsolve
        from scipy.linalg import norm
        from scipy.integrate import solve_ivp
        from matplotlib.pyplot import figure

        if n0 is None:
            n0 = self.n0()
        if ext is None:
            ext = zeros((len(n0),))
        #Te /= 1.602e-19
        #ne /= 1e6
        # Set volume to be cm**3
        #vol *= 1e6
        # Set initial density in [1/cm**3]
        # Get the rate matrices 
        mat, extm = self.M(Te, ne, Ti, ni, E, write=True) 
        ext += extm
        if gl is True:
            # Get matrices and eigenvalues
            mat, ext, n = self.gl_crm(mat, ext, n=n) 
        # Set external sources to be [1/cm**3/s]
#        ext[0] = (psorgc + diva + srca + bga)/vol
#        ext[1] = (divm + srcm + bgm)/vol
        # Normalize ionzation rate to UEDGE values 
#        mat[0,0] = (psor/ne)/vol
        mat[0,0] += ioniz
        if plot is True:
            # Simulate to SS
            ss_sol = solve_ivp(lambda x,y: self.ddt(y, x, mat, ext), (0, 100), 
                    n0, 'LSODA', dense_output=True)
            f = figure(figsize = (7,7))
            ax = f.add_subplot(111)
            for i in range(len(n0)):
                line = '-' + '-'*('H2' in self.slist[i])
                ax.semilogx(ss_sol.t, ss_sol.y[i,:], line, label=self.slist[i])
                if i == 0:
                    ax.semilogx(ss_sol.t, ss_sol.y[i,:], 'k.', 
                            label=self.slist[i]) 
            ax.legend(ncol=3)
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Density [cm**-3]')
            f.show()
        if dt is False:
            return fsolve(self.ddt, n0, args=(0, mat, ext))  
        else: 
            return solve_ivp(lambda x,y: self.ddt(y, x, mat, ext), (0, 100), 
                    n0, 'LSODA', dense_output=True).y[:,-1]


    def gl_crm(self, mat, ext, Sext=True, n0=None):
        ''' Returns the P-space matrices according to Greenland 2001

        Parameters
        ----------
        mat : ndarray
            2D rate matrix
        ext : ndarray
            1D external sink/source vector
        Sext : boolean, optional (default: True)
            switch for considering 'external' sinks/sources (e.g. 
            re-association)
        n0 : list of floats, optional (default: None)
            Initial density distribution of the species defined in the 
            input. List must be of the same length as the species list,
            and the order is determined by the input file. Defaults to 
            the initial densities defined in the input if n0 is None.

        # TODO: Eigenvalues into its own function!
        Returns 
        -------
        Meff, GPp, nP0p if eigen is False
        Meff : ndarray
            the Greenland Np x Np effective rate matrix for the 
            P-species
        GPp : ndarray
            the Greenland Np effective external sink/source vector for
            the P-species
        nP0p : ndarray
            the Greenland Np effective initial densities for the 
            P-species
        eigenvectors, eigenvalues if eigen is True
        eigenvectors : ndarray
            The eigenvectors of the supplied rate matrix
        eigenvalues : ndarray
            The eigenvalues of the supplied rate matrix
        '''
        from numpy import matmul, diag, real
        from numpy.linalg import inv, eig

        # Use input n0 as default unless explict n0 requested
        if n0 is None: 
            n0 = self.n0() 
        # Create block matric from M
        MP = mat[:self.Np, :self.Np]
        MQ = mat[self.Np:, self.Np:]
        V = mat[self.Np:, :self.Np]
        H = mat[:self.Np, self.Np:]
        # Calculate Meff
        Meff = (MP - matmul(matmul(H, inv(MQ)), V))
        # Diagonalize M
        eigs, T = eig(mat)
        # Order the eigenvalues and vectors in increasing magnitude
        # TODO should this be done for all the eigenvectors or separately   
        # between the Q and P space??
        # --> Does not seem to affect the solution
        eigind = abs(eigs[self.Np:]).argsort()[::1] 
        e1 = eigs[self.Np:]
        T1 = T[:,self.Np:]
        e1 = e1[eigind]
        T1 = T1[:,eigind]
        #eigs=eigs[eigind]
        #T=T[:,eigind]
        eigs[self.Np:] = e1
        T[:,self.Np:] = T1
        # Only use real parts TODO is this OK??
        eigs = real(eigs)
        D = diag(eigs)
        # Create block matrix
        TQ = T[self.Np:, self.Np:]
        TP = T[:self.Np, :self.Np]
        delta = T[self.Np:, :self.Np]
        Delta = T[:self.Np, self.Np:]
        # Calculate P-space CRM
        GPp = real((Sext is True)*ext[:self.Np] - matmul(matmul(Delta,
                inv(TQ)), ext[self.Np:]))
        # TODO: Use n0 instead of nP0P?
        nP0p = real(n0[:self.Np] - matmul(matmul(Delta, inv(TQ)), 
                n0[self.Np:])) 
        #nP0p = n0[:self.Np]
        return Meff, GPp, nP0p

  
    def ddt(self, n, t, mat, ext):
        ''' RHS of the linearized system of ODEs ''' 
        from numpy import matmul

        return matmul(mat, n) + ext

    def evaluate_CRM(self, Te, ne, Ti=None, ni=None, E=0.1, printout=True):
        ''' **INCOMPLETE** Evaluates the current CRM  according to Greenland 2001

        Evaluates the separation between the minimum lifetimes and CRM
        parameters of the P-space and that of the Q-space, which need
        to be sufficiently far apart for the CRM approximation to be 
        accurate.

        Parameters
        ----------
         Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        Ti : float, optional (default: None)
            ion temperature in eV. Default sets Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Default sets ni=ne
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        printout : boolean, optional (default: True)
            switch to print results (True) or return the results (False)

        Returns
        -------
        None if printout is True
        maxdelta, maxnorm, tauPmin, tauQmax
        maxdelta : float
            maximum norm of the Greenland delta-parameter
        maxnorm : float
            maximum of the norm of inv(Q_t)*delta
        tauPmin : float
            minimum lifetime of the P-species
        tauQmax : float
            maximum lifetome of the Q-species
        '''
        from numpy import matmul
        from numpy.linalg import inv

        # Get matrices and eigenvalues/-vectors
        M, T, D = self.gl_crm(*self.M(Te, ne, Ti, ni, E, write=False), 
                matrices=True) 
        TQ, delta = T[self.Np:, self.Np:], T[self.Np:, :self.Np]
        tau = 1/abs(D)
        tauPmin = min(tau[:self.Np])
        tauQmax = max(tau[self.Np:])
        norm = matmul(inv(TQ), delta)
        maxnorm = max([sum(abs(norm[:,i])) for i in range(norm.shape[1])])
        maxdelta = max([sum(abs(delta[:,i])) for i in range(delta.shape[1])])
        if printout is True:
            print("Validity of CRM with Np={}".format(self.slist[:self.Np]))
            print("    max(|delta|)={:.2E}".format(maxdelta))
            print("    max(|inv(T_Q)*delta|)={:.2E}".format(maxnorm))
            print("    min(tau_P)={:.2E}".format(tauPmin))
            print("    max(tau_Q)={:.2E}".format(tauQmax))
        else:
            return [maxdelta ,maxnorm, tauPmin, tauQmax]

    
    def generate_CRM(self, Te, ne, kappa, Ti=None, ni=None, E=0.1, epsilon=1):
        ''' **INCOMPLETE** Generates the optimal CRMs per Greenland 2001 

        Algorithm for finding the optimal set of P-species for a given 
        temperature and density. 
        
        Parameters
        ----------
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        kappa : float
            the threshold of the algorithm, <<1
        Ti : float, optional (default: None)
            ion temperature in eV. Default sets Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Default sets ni=ne
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        epsilon : float, optional (default: 1)
            sensitivity??

        Returns
        -------
        None
        '''
        from numpy import array,where,matmul
        from numpy.linalg import inv

        # Get matrices and eigenvalues/-vectors
        M, T, D = self.gl_crm(*self.M(Te, ne, Ti, ni, E, write=False),
                matrices=True) 
        # Construct indicator matrix
        I = abs(T)
        I[I <= kappa] = 0
        I[I > kappa] = 1
        # TODO: What is epsilon?
        # Create Z-array
        Z = []
        for i in I: # Loop over species
            N = 0
            while i[N] == 0: # Count leading zeros
                N += 1
                # Break if all species are 0
                if N == len(i):
                    break
            Z.append(N)
        # Get increasing order of leading zeros
        Z = array(Z)
        Zind = Z.argsort()[::1] 
        # Check whether indicator matrix is singular
        if max(Z) == M.shape[1]:
            print('WARNING! Matrix I is singular')
        # Identify a boolean array for valid Np
        Np = [Z[Zind][i] >= i for i in range(len(Z))]
        di = range(len(Z)) - Zind # Track change in index space
        sorted_D = D[Zind]
        sorted_species = [self.slist[i] for i in Zind]
        rej = []
        for i in range(len(Z)):
            if Np[i] == True:
                if i + 1 == len(Z): 
                    continue # Full model, uninteresting
                TQ = T[i+1:,i+1:]
                delta = T[i+1:,:i+1]
                norm = matmul(inv(TQ), delta)
                maxnorm = max([sum(abs(norm[:,i])) for i in 
                        range(norm.shape[1])])
                if maxnorm < epsilon:
                    print('Possible CRM NP={} model with TQinv*delta={:.2E}, '
                            'tauP={:.2E} s, tauQ={:.2E} s, and  P-space {} '
                            'found.'.format(i + 1, maxnorm,
                            1/abs(min(sorted_D[:i + 1])),
                            1/abs(max(sorted_D[i + 1:])), 
                            sorted_species[:i + 1])) 
                else:
                    rej.append(maxnorm)
        if len(rej) > 0:
            print('Minimum rejected maxnorm={:.2E}'.format(min(rej)))

        
    def intensity(
            self, Te, ne, vol, E=0.1, Ti=None, ni=None, Sext=True, n0=None, 
            ext=None, ioniz=0, units='l',write=False, norm=False):
        ''' Returns energy and gammas and the respective counts at SS

        Parameters
        ----------
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        vol : float
            volume of the domain in cm**3
        Ti : float, optional (default: None)
            ion temperature in eV. Default sets Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Default sets ni=ne
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        Sext : boolean, optional (default: True)
            switch for considering 'external' sinks/sources (e.g. 
            re-association)
        n0 : ndarray, optional (default: None)
            Initial density distribution of the species defined in the 
            input. List must be of the same length as the species list,
            and the order is determined by the input file. Defaults to 
            the initial densities defined in the input if n0 is None.
        ext : ndarray, optional (default: None)
            External sinks and sources for the species, not considered 
            in CRUMPET, including transport, artificial sinks and 
            sources, etc. Must be of length Nx1, as defined by the input
            species. Defaults to zeros, given in cm**-3 s**-1
        ioniz : float, optional (default: None)
            Additional ionization source (<0) to applied to the first 
            species. Necessary for evaluating CRM-transport code 
            couplings where ionization is handled by the auxilliary
            code, and not included in the CRUMPET model. Given in 
            cm**-3 s**-1
        units : string, optional (default: 'l')
            units of the plot ordinate
            ev - power in eV/s
            v - wavenumber cm**-1
            f - frequency Hz
            l - wavelength nm
            Å - wavelength Å
        write : boolean, optional (default: False)
            write intensity matrices to logs
        norm : boolean, optional (default: False)
            normalize the lines to the total as determined by units

        Returns
        -------
        intensities : len-2 tuple of ndarray
            the first ndarray contains atomic line intensities, and the
            second molecular line intensities. The first entry along
            axis=0 is are the gamma energies in the units requested, and 
            the second are counts in s**-1
        ''' 
        from numpy import reshape, zeros, multiply, where
        from numpy import transpose, matmul, array, mean, all, diag
        from scipy.optimize import minimize
        from scipy.integrate import solve_ivp
        from matplotlib.pyplot import figure


        n_ss = self.steady_state(Te, ne, vol, E, Ti, ni, Sext, False, False, 
                n0, ext, ioniz)
        x = self.E(Te, ne, Ti, Te, E, write=False)
        y = self.I(Te, ne, Ti, Te, E, write=False)
        if write:
            self.write_matrix(y[0][0], n, 'Ia', Te, ne, Ti, ni, E, 
                    form='{:1.2E}')
            self.write_matrix(y[1][0], n, 'Im', Te, ne, Ti, ni, E,
                    form='{:1.2E}')
        ret = []
        for i in range(2):
            arr = zeros((len(n_ss)**2, 2))
            # Use the steady-state solution to find transitions 
            Y = matmul(y[i][0], diag(n_ss))
            arr[:,0] = reshape(x[i][0], (1, len(n_ss)**2))
            arr[:,1] = reshape(Y, (1, len(n_ss)**2))
            # Drop non-radiative entries
            arr = transpose(arr[~all(abs(arr) < 1e-8, axis=1)])
            if units == 'ev':     
                arr[0,:] = arr[0,:]
            elif units == 'v':    
                arr[0,:]=1e7/(1239.84193/arr[0,:])
            elif units == 'f':
                arr[0,:]=241.798*arr[0,:]
            elif units == 'l':
                arr[0,:]=1239.84193/arr[0,:]
            elif units == 'Å':
                arr[0,:]=12398.4193/arr[0,:]
            if norm is True:
                nrm = sum(arr[1,:])
            else:
                nrm = 1
            arr[1,:] = arr[1,:]/nrm
            ret.append(arr)
        return ret
            

    def species_rate(self,species,Te,ne,Ti=None,ni=None,E=0.1):
        ''' Returns the total reaction  rate of a given species

        Returns the net reaction rate of all reactions that have
        the specified species handle as a reactant. For CRM species,
        this is analog to depletion rate, but for background species
        it is net rate of reactions involving the background species.

        Parameters
        ----------
        species : string
            species handle as defined in the input file
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        Ti : float, optional (default: None)
            ion temperature in eV. Default sets Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Default sets ni=ne
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        
        Returns
        -------
        totrate : float
            the net rate of all reactions having species as a reactant 
        '''

        if Ti is None: 
            Ti = Te
        if ni is None: 
            ni = ne
        ret = 0
        for dkey, database in self.reactions.items():
            mult = 1
            for rkey, r in database.items():
                if species in r.reactants:
                    if r.type == 'COEFFICIENT':
                        mult = 1/ne
                    ret += r.rate(Te, Ti, E=E, ne=ne)*mult
        return ret


    def n0(self):
        ''' Returns initial densities as specieifed in the input '''
        from numpy import zeros
        ret = zeros((len(self.species),))
        for s in range(len(self.species)):
            try:
                ret[s] = self.species[self.slist[s]]['n'] 
            except:
                pass
        return ret


    def get_reaction(self, database, name):
        ''' Returns the Reaction object for the requested reaction
        
        Parameters
        ----------
        database : string
            the database where to look for the reaction
        name : string
            the reaction name/handle to look for
        '''
         
        return self.reactions[database][name]


    def get_rate(self, database, name, T, n, E=0.1):
        ''' **INCOMPLETE** Returns the rate of the requested reaction '''
        return self.get_reaction(database,name).rate(T, T, E, n)


    def totparticles(self, arr):
        ''' Returns the total number of nuclei for a vector or matrix

        Parameters
        ----------
        arr : ndarray
            array or matrix of species densities along axes=0

        Returns
        -------
        totnuclei : float
            total number of nulclei in the system
        '''
        from numpy import zeros

        try:    # If matrix 
            ret = zeros((arr.shape[1],))
            for j in range(arr.shape[1]):   # Loop through eac time-step
                for i in range(len(arr)):   # Loop through each species in arr
                    ret[j] += (1 + ('{}2'.format(self.isotope) 
                            in self.slist[i]))*arr[i,j]
        except: # If vector
            ret = 0
            for i in range(len(arr)): # Loop through each species in arr
                ret += (1 + ('{}2'.format(self.isotope) 
                        in self.slist[i]))*arr[i]
        return ret

    def totmol(self, arr):
        from numpy import zeros
        ''' Returns the total number of molecules in a vector or matrix

        Parameters
        ----------
        arr : ndarray
            array or matrix of species densities along axes=0

        Returns
        -------
        totnuclei : float
            total number of molecules in the system
        '''
        try:    # If matrix 
            ret = zeros((arr.shape[1],))
            for j in range(arr.shape[1]):   # Loop through eac time-step
                for i in range(len(arr)):   # Loop through each species in arr
                    ret[j] += (('{}2'.format(self.isotope) 
                            in self.slist[i]))*arr[i,j]
        except: # If vector
            ret = 0
            for i in range(len(arr)): # Loop through each species in arr
                ret += (('{}2'.format(self.isotope) 
                        in self.slist[i]))*arr[i]
        return ret

    def get_species_ind(self,species):
        ''' TODO: documentation'''
        return list(self.species).index(species)














    def populate_old(self,mode,Te,ne,Ti=None,ni=None,E=0,Sind=None,Tm=0,Iind=0):
        ''' **REDUNDANT*** '''
        from numpy import zeros,array,sum,transpose
        

        if mode=='diagnostic':
            # Setup a 2D diagnostic list for the matrix and a list for 
            # the external source
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


    def S(self, Te, ne, Ti=None, ni=None, E=0.1,
                                Tm=0,write=False):
        ''' **OBSOLETE** '''
        if Ti is None: Ti=Te # Check for Ti, set if necessary
        if ni is None: ni=ne # Check for ni, set if necessary

        mat,ext=self.populate('Sgl',Te,ne,Ti,ni,E,Tm=Tm)    


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




    def full_nt(self,Te,ne,t,gl,Ti=None,ni=None,E=0.1,n=None,Sext=True):
        ''' **OBSOLETE** '''
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
        return solve_ivp(lambda x,y: self.solve_crm(x,y,mat,ext),
                                            (0,t),n,method='LSODA') 



        
    def dEdt_old(self,t,Te,ne,Ti=None,ni=None,E=0.1,Tm=0,
                Sext=True,write=False,gl=True,n=None,Qres=True):
        ''' **OBSOLETE** '''
        from numpy import block,zeros,matmul,reshape,sum
        from  numpy.linalg import inv
        from scipy.integrate import solve_ivp


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
        
            U=self.Sgl(Te,ne,Ti,ni,E,Tm,write=write)


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
                U=self.S(Te,ne,Ti,ni,E,Tm,write=True)
    

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
                U=self.S(Te,ne,Ti,ni,E,Tm,write=True)

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

        return solve_ivp(lambda x,y: self.solve_crm(x,y,mat,ext),(0,t),
                                            n,'LSODA',dense_output=True)

        


