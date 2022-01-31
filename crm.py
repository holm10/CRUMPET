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


    def __init__(self, species, bg, reactions, path='.', 
                rdata=None, settings = {}, writelog=False):
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
            'educts' : list of strings
                each element is a handle of the educts. At least one
                of the species must be a background species, and only
                two educts are currently supported
            'products' : list of strings
                each element is a handle of the products
            'K' : string
                a string of the net kinetic energy transfer for the 
                background educts to the products. Presently, all
                energy is assumed to end up as ion/atom kinetic energy.
                Supports calculations in the string, e.g. '2*3+4'
        verbose : boolean
            switch whether to print the CRM setup to the prompt: the 
            same info is written to {Crm.path}/logs regardless
        NP : int (default: 1)
            the number of P-species to be considered in the model: the 
            NP first entries in species (e.g. the input file) are taken
            as the P-species. The remainder are assumed Q-species.
        path : string, optional (default: '.')
            a path to which all subsequent paths are relative to. 
            Relative or an absolute path.
        vmax : int, optional (default: 14)
            number of vibrational excited levels to be considered for 
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
        from numpy import ones, zeros
        from CRUMPET.reactions import Reaction
        from copy import deepcopy

        # Store class objects
        self.species = species
        self.slist = list(self.species)
        self.N = len(self.slist)
        self.bg = bg
        self.path = path
        self.reactions = {} 

        # Optional settings 
        # Set defaults
        self.isotope = 'H'
        self.mass = 1
        self.verbose = False
        self.Np = 1
        self.vmax = 14
        self.nmax = 8
        self.source = zeros((self.N,))
        for s in self.slist:
            self.species[s]['n'] = 0
        self.evolvearr = ones((self.N,))
        self.evolvemat = ones((self.N, self.N))

        for key, value in settings.items():
            if key.upper() == 'ISOTOPE':
                self.isotope = value
            elif key.upper() == 'MASS':
                self.mass = int(value)
            elif key.upper() == 'VERBOSE':
                self.verbose = (value.upper() not in ['0', 'FALSE', 
                                'NO', 'OFF'])
            elif key.upper() == 'NP':
                self.Np = int(value)
            elif key.upper() == 'VMAX':
                self.vmax = int(value)
            elif key.upper() == 'NMAX':
                self.nmax = int(value)
            elif key.upper() == 'N0':
                for s in value:
                    self.species[s.split()[0]]['n'] = float(s.split()[1])
            elif key.upper() == 'STATIC':
                for s in value:
                    self.evolvearr[self.slist.index(s)] = 0
            elif key.upper() == 'SOURCE':
                for s in value:
                    [spec, val] = s.split()
                    self.source[self.slist.index(spec)] = float(val)
            else:
                print('Setting "{}" not recognized! Ignoring...'.format(key))
        
        for i in range(self.N):
            self.evolvemat[i] *= self.evolvearr[i]

        '''
        for database, subdatabase in reactionlist.items():
            for h123, reactions in subdatabase.items():
                for rname, data in reactions.items():
                    print(database, h123, rname, data)
        '''
        def set_coeffs(db, h123, r, rdata):
            if db in ['AMJUEL','HYDHEL','H2VIBR']:
                try:
                    return rdata[db][h123][r.upper()]
                except:
                    if '-' in r:
                        return rdata[db][h123][r.split('-')[0].upper()]
                    else: 
                        print('Reaction "{}" of type "{}" not '
                            'found in database "{}"!'.format(r, h123, db))
                    return False
            else:
                return False
            
        
        self.reactions = deepcopy(reactions)
        for db, dbcontent in reactions.items():
            if (db in ['AMJUEL', 'HYDHEL', 'H2VIBR']) and (db not in rdata.keys()):
                print('Required database "{}" not found, aborting!'.format(db))
                return
            for h123, h123content in dbcontent.items():
                for rea, reactiondata in h123content.items():
                    # Split the definitions into names and databases
                    reastr = reactiondata[0].upper()
                    if db in ['FCF', 'AIK']:
                        ''' Branch for vibrational H2 transitions by factors '''
                        del self.reactions[db][h123][rea]
                        datafile = reactiondata[1]
                        if db == 'AIK':
                            [ve, vp] = [int(x) for x in reactiondata[2].split()]
                            # Read table of Einstein coefficients
                            for x in range(ve+1):
                                for y in range(vp+1):
                                    reabuff = self.XY2num(rea,x,y)
#                                    print(rdata[db][datafile][y,x])
                                    self.reactions[db][h123][reabuff] = \
                                            Reaction(db, h123, 
                                            reabuff, [self.XY2num(reactiondata[0],x,y)], 
                                            rdata[db][datafile][y,x],
                                            self.bg, self.species, 
                                            self.isotope, self.mass)
                        elif db == 'FCF':
                            # Read table of Franck-Condon factors
                            vl = [int(x) for x in reactiondata[2].split()]
                            vu = [int(x) for x in reactiondata[3].split()]
                            if len(vl) == 1:
                                x = vl[0]
                                # Only upper state to be read from FCFs
                                for y in range(vu[0], vu[1]+1):
                                    reabuff = self.XY2num(rea,y)
                                    self.reactions[db][h123][reabuff] = \
                                            Reaction(db, h123, 
                                            reabuff, [self.XY2num(reactiondata[0],y)]+reactiondata[4:], 
                                            0,
                                            self.bg, self.species, 
                                            self.isotope, self.mass, fcf=rdata[db][datafile][x,y])
                            elif len(vu) == 1:
                                # Only lower state to be read from FCFs
                                y = vu[0]
                                for x in range(vl[0], vl[1]+1):
                                    reabuff = self.XY2num(rea,y)
                                    self.reactions[db][h123][reabuff] = \
                                            Reaction(db, h123, 
                                            reabuff, [self.XY2num(reactiondata[0],y)]+reactiondata[4:], 
                                            0,
                                            self.bg, self.species, 
                                            self.isotope, self.mass, fcf=rdata[db][datafile][x,y])
                                # Only lower state to be read from FCFs
                            else:
                                # Both upper and lower states are to be read
                                for x in range(vl[0], vl[1]+1):
                                    for y in range(vu[0], vu[1]+1):
                                        reabuff = self.XY2num(rea,x,y)
                                        self.reactions[db][h123][reabuff] = \
                                                Reaction(db, h123, 
                                                reabuff, [self.XY2num(reactiondata[0],x,y)]+reactiondata[4:], 
                                                0,
                                                self.bg, self.species, 
                                                self.isotope, self.mass, fcf=rdata[db][datafile][x,y])
    
        
                    else:
                        ''' Branch for ground-state transitions '''
                        # Loop through states, if necessary. Dynamicallt set boundaries
                        # according to electronic or vibrational transitions
                        isn, isv  = ('N=$' in reastr), ('V=$' in reastr)
                        if isn or isv:
                            del self.reactions[db][h123][rea]
                        for x in range(isn, 1 + isv*self.vmax + isn*self.nmax):
                            # Vibrational/electronic dependence present
                            if '$' in reastr:
                                buff = reactiondata.copy()
                                reabuff = rea
                                # Substitute state into educts and product strings 
                                buff[0] = self.XY2num(buff[0], x)
                                reabuff = self.XY2num(reabuff, x)
                                # %%% n- or v-dependent upper rate %%%
                                if '&' in reastr:
                                    # Relaxation - only move down during reaction
                                    if h123.upper() == 'RELAXATION':
                                        for y in range(int(isn), x):
                                            rebuff = buff.copy()
                                            rereabuff = reabuff
                                            rebuff[0] = self.XY2num(rebuff[0], x, y) 
                                            rereabuff = self.XY2num(rereabuff, x, y) 
                                            coeffs = set_coeffs(db, h123, rereabuff.upper(), rdata)
                                            if coeffs is False:
                                                coeffs = (x, y)
                                            self.reactions[db][h123][rereabuff] = \
                                                    Reaction(db, h123, 
                                                    rereabuff, rebuff, coeffs,
                                                    self.bg, self.species, 
                                                    self.isotope, self.mass)
                                    # Excitation - only move up during reaction
                                    elif h123.upper() == 'EXCITATION':
                                        for y in range(x + 1, 1 + isn*self.nmax + 
                                                            isv*self.vmax):
                                            rebuff = buff.copy()
                                            rereabuff = reabuff
                                            rebuff[0] = self.XY2num(rebuff[0], x, y) 
                                            rereabuff = self.XY2num(rereabuff, x, y) 
                                            coeffs = set_coeffs(db, h123, rereabuff.upper(), rdata)
                                            if coeffs is False:
                                                coeffs = (x, y)
                                            self.reactions[db][h123][rereabuff] = \
                                                    Reaction(db, h123, 
                                                    rereabuff, rebuff, coeffs,
                                                    self.bg, self.species, 
                                                    self.isotope, self.mass)
                                    # Ladder-like process assumed
                                    elif 'H.' in h123.upper():
                                        for y in range(-1*('&' in buff[0]), 2, 2):
                                            # Limit transitions to [0,vmax]
                                            if x + y in range(self.vmax + 1):
                                                rebuff = buff.copy()
                                                rereabuff = reabuff
                                                # Retain intial and final states in name
                                                rebuff[0] = self.XY2num(rebuff[0], x, x + y) 
                                                rereabuff = self.XY2num(rereabuff, x, x + y) 
                                                coeffs = set_coeffs(db, h123, rereabuff.upper(), rdata)
                                                if coeffs is False:
                                                    coeffs = (x, y)
                                                self.reactions[db][h123][rereabuff] = \
                                                        Reaction(db, h123, 
                                                        rereabuff, rebuff, coeffs,
                                                        self.bg, self.species, 
                                                        self.isotope, self.mass)
                                    else:
                                        print('Reaction specifier "{}" not' 
                                            ' recognized. Omitting '
                                            '{} {} {}'.format(h123, db, h123, rea))
                                # Independent upprt state
                                else:
                                    coeffs = set_coeffs(db, h123, reabuff.upper(), rdata)
                                    if coeffs is False:
                                        coeffs = x
                                    self.reactions[db][h123][reabuff] = \
                                            Reaction(db, h123, reabuff, buff,
                                            coeffs, self.bg, self.species, 
                                            self.isotope, self.mass)
                            else:
#                                if db.upper() == 'USER':
 #                                   if h123.upper() == 'COEFFICIENT':
  #                                      print(reactiondata, rea)
   #                                     self.reactions[db][h123][rea] = Reaction(db, h123, rea,
    #                                                reactiondata, coeffs, self.bg, 
     #                                               self.species, self.isotope, self.mass)
      #                          else:
                                coeffs = set_coeffs(db, h123, rea.upper(), rdata)
                                if coeffs is False:
                                    coeffs = None
                                self.reactions[db][h123][rea] = Reaction(db, h123, rea,
                                            reactiondata, coeffs, self.bg, 
                                            self.species, self.isotope, self.mass)


        # Define reactions for UEDGE raditation
        #self.ionizrad=Reaction('IONIZRAD','UE',
                #self.get_coeff('UE','IONIZRAD'),'UE',['','',''],bg,species,
                #['','',None,None,[0,0,0,0]])
        #self.recrad=Reaction('RECRAD','UE',self.get_coeff('UE','RECRAD'),
                #'UE',['','',''],bg,species,['','',None,None,[0,0,0,0]])
        #return
        # Create the rate matrix M based on the input
        print('Constructing functional rate matrix')
        self.create_M()
        print('Constructing functional emissivity matrix')
        self.create_EI()
        print('Constructing functional energy transfer matrix')
        self.create_Sgl()


        # Ensure that there is a logs directory under the run path
        try:
            mkdir('{}/logs'.format(self.path))
        except:
            pass
        if writelog is True:
            # Write a log of the CRM setup path/logs
            with open('{}/logs/setup.log'.format(self.path), 'w') as f:
                f.write('CRUMPET run in {} on {}\n'.format(getcwd(),
                        str(datetime.now())[:-7]))
                f.write('Defined species:\n')
                for i in self.slist:
                    f.write('    {}\n'.format(i))
                f.write('Defined reactions:\n')
                for dkey, database in self.reactions.items():
                    for hkey, h123 in database.items():
                        for rkey, reaction in h123.items():
                            f.write('{}\n'.format(reaction.print_reaction()))
        # Output to stdout if run verbosely
        if self.verbose:
            with open('logs/setup.log', 'rt') as f:
                for m in f:
                    print(m.strip())
        # Do the same for a Diagnostic rate matrix displaying 
        # reaction correlations
        self.DIAGNOSTIC()


    def create_EI(self):
        ''' Function populating and returning the rate matrix and source

        '''
        from numpy import zeros,array,sum,transpose
        from math import isinf
        from tqdm import tqdm
        
        N = len(self.species)
        E = []
        I = []
        GE = []
        GI = []
        for i in range(N):
            E.append([])
            I.append([])
            GE.append([])
            GI.append([])
            for j in range(N):
                E[i].append([])
                I[i].append([])

        for i in tqdm(range(N)):
            #%%% Walk through each row (species) %%%
            for dkey, database in self.reactions.items():
                for h123, h123data in database.items():
                    for rkey, r in h123data.items():
                        #%%% Sort the species of each reaction into %%%
                        #%%%  the appropriate column %%%
                        for frag in range(len(r.products)):    
                            ''' SOURCE '''
                            # If fragment enters row, add in appropriate column as
                            # defined by reactant
                            if (r.products[frag] in self.slist) and (
                                    self.slist.index(r.products[frag]) == i):
                                # External flag triggered, store to external source
                                common = set(r.educts).intersection(\
                                        set(self.slist))
                                if len(common) == 0:
                                    ''' EXTERNAL SOURCE '''
                                    # TODO: what reactions end up here? External contributions to radiation???
                                    #if (r.e is not False) and (r.p is not False):
                                    GI[i].append(r.rate)
                                    GE[i].append(r.getS)
                                # No trigger of external source, store to 
                                # appropriate location in matrix
                                else: 
                                    ''' INTERNAL SOURCE '''
                                    j = self.slist.index(common.pop())
                                    # TODO: Other check than positive energy change?
                                    # How should "spontaneous excitation" be considered?
                                    I[i][j].append(r.rate)
                                    E[i][j].append(r.getS)
        self.I = I
        self.E = E
        self.GI = GI
        self.GE = GE



    def eval_EI(self, Te, ne, Ti, ni, Energy, Tm):
        ''' Returns the reaction rate matrix for the local plasma parameters '''
        from numpy import zeros
        E, GE = zeros((self.N, self.N, 2)), zeros((self.N, 2))
        I, GI = zeros((self.N, self.N, 2)), zeros((self.N, 2))
        for i in range(self.N):
            for j in range(self.N):
                for r in range(len(self.E[i][j])):
                    S = self.E[i][j][r](Te, Ti, Tm, Energy, marker=True)[-2:,0]
                    E[i,j] += S*self.E[i][j][r](Te,Ti,Tm,Energy,marker=False)[-2:,0]
                    I[i,j] += S*self.I[i][j][r](Te=Te, ne=ne, Ti=Ti, ni=ni, E=Energy, frequency=True)
            for r in range(len(self.GE[i])):
                S = self.GE[i][r](Te, Ti, Tm, E)[-2:,0]
                GE[i] += S* self.GE[i][r](Te, Ti, Tm, E, marker=False)[-2:,0]
                GI[i] += S*self.GI[i][r](Te=Te, ne=ne, Ti=Ti, ni=ni, E=Energy, frequency=True)
        
        return E, GE, I, GI




    def create_Sgl(self):
        ''' Function populating and returning the rate matrix and source

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
        Tm : float, optional (default: None)
            local molecular temperature for calculating energy losses 
            associated with molecules. Defaults to 2E/3 if Tm=None.

        Returns
        -------
        mat : ndarray
            Rate matrix containing the rates at the given 
            temperature and density. A 2D list in case mode=diagnostic
        rec : ndarray
            Vector of recombination sources. List if mode=diagnostic
        '''
        from numpy import zeros,array,sum,transpose
        from math import isinf
        from tqdm import tqdm
        
        N = len(self.species)
        S = []
        SS = []
        Sfac = []
        GS = []
        GSS =[]
        for i in range(N):
            S.append([])
            Sfac.append([])
            SS.append([])
            GS.append([])
            GSS.append([])
            for j in range(N):
                S[i].append([])
                Sfac[i].append([])
                SS[i].append([])

        count=0
        for i in tqdm(range(N)):
            #%%% Walk through each row (species) %%%
            for dkey, database in self.reactions.items():
                for h123, h123data in database.items():
                    for rkey, r in h123data.items():
                        #%%% Sort the species of each reaction into %%%
                        #%%%  the appropriate column %%%
                        # Loop through the reaction products
                        for frag in range(len(r.products)):    
                            ''' SOURCE '''
                            # Do nothing if background fragment
                         #   if r.products[frag] not in self.slist: 
                         #       continue
                            # If fragment enters row, add in appropriate column as
                            # defined by reactant
                            if (r.products[frag] in self.slist) and (
                                    self.slist.index(r.products[frag]) == i):
                                # External flag triggered, store to external source
                                common = set(r.educts).intersection(\
                                        set(self.slist))
                                if len(common) == 0:
                                    ''' EXTERNAL SOURCE '''
                                    GS[i].append(r.rate)
                                    GSS[i].append(r.getS)
                                # No trigger of external source, store to 
                                # appropriate location in matrix
                                else: 
                                    ''' INTERNAL SOURCE '''
                                    count +=1
                                    ind = self.slist.index(common.pop())
                                    S[i][ind].append(r.rate)
                                    Sfac[i][ind].append(r.get_nsum)
                                    SS[i][ind].append(r.getS)
        self.S = S
        self.Sfac = Sfac
        self.SS = SS
        self.GS = GS
        self.GSS = GSS

    def eval_Sgl(self, Te, ne, Ti, ni, Energy, Tm):
        ''' Returns the reaction rate matrix for the local plasma parameters '''
        from numpy import zeros
        S, GS = zeros((self.N, self.N, 5)), zeros((self.N, 5))
        for i in range(self.N):
            for j in range(self.N):
                for r in range(len(self.S[i][j])):
                    Ssource = self.SS[i][j][r](Te, Ti, Tm, Energy)
                    S[i,j] += self.Sfac[i][j][r](ne,ni)*self.S[i][j][r](Te=Te, ne=ne, Ti=Ti, ni=ni, E=Energy, frequency=False)*Ssource[:,0] + Ssource[:,1]
            for r in range(len(self.GS[i])):
                Ssource = self.GSS[i][r](Te, Ti, Tm, Energy)
                GS[i] += self.GS[i][r](Te=Te, ne=ne, Ti=Ti, ni=ni, E=Energy, frequency=True)*Ssource[:,0] + Ssource[:,1]
        
        return S, GS


    """
    def eval_M_para(self, Mcol, Te, ne, Ti, ni, E):
        ''' Returns the reaction rate matrix for the local plasma parameters '''
        from numpy import zeros
        M = zeros((self.N))
        for j in range(self.N):
            for r in range(len(Mcol[j])):
                M[j] += Mcol[j][r](Te=Te, ne=ne, Ti=Ti, ni=ni, E=E, frequency=True)
        return M

    def M_para(self, Te, ne, Ti, ni, E):
        from joblib import Parallel, delayed
        from multiprocessing import cpu_count
        
        ret = Parallel(n_jobs=cpu_count())(delayed(self.eval_M_para)(Mcol, Te, ne, Ti, ni, E) for Mcol in self.M)

        return ret
    """

    def eval_M(self, Te, ne, Ti, ni, E):
        ''' Returns the reaction rate matrix for the local plasma parameters '''
        from numpy import zeros
        M, G = zeros((self.N, self.N)), zeros((self.N))
        for i in range(self.N):
            for j in range(self.N):
                for r in range(len(self.M[i][j])):
                    M[i,j] += self.M_mult[i][j][r]*self.M[i][j][r](Te=Te, ne=ne, Ti=Ti, ni=ni, E=E, frequency=True)
            for r in range(len(self.G[i])):
                G[i] += self.G_mult[i][r]*self.G[i][r](Te=Te, ne=ne, Ti=Ti, ni=ni, E=E, frequency=True)
        return M, G


    def eval_R(self, Te, ne, Ti, ni, E):
        ''' Returns the reaction rate matrixes NOT weighted by density'''
        from numpy import zeros
        # NOTE: How to treat radiative decays? Different units!
        R, GR = zeros((self.N, self.N)), zeros((self.N))
        for i in range(self.N):
            for j in range(self.N):
                for r in range(len(self.M[i][j])):
                    R[i,j] += self.M_mult[i][j][r]*self.M[i][j][r](Te=Te, ne=ne, Ti=Ti, ni=ni, E=E, frequency=False)
            for r in range(len(self.G[i])):
                GR[i] += self.G_mult[i][r]*self.G[i][r](Te=Te, ne=ne, Ti=Ti, ni=ni, E=E, frequency=False)
        return R, GR

    def create_M(self):
        ''' Function populating and returning the rate matrix and source

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

        Returns
        -------
        mat : ndarray
            Rate matrix containing the rates at the given 
            temperature and density. A 2D list in case mode=diagnostic
        rec : ndarray
            Vector of recombination sources. List if mode=diagnostic
        '''
        from numpy import zeros,array,sum,transpose
        from math import isinf
        from tqdm import tqdm
        
        N = len(self.species)
        # Setup a matrix and vector
        M = []
        G = []
        G_mult = []
        M_mult = []
        for i in range(N):
            G.append([])
            G_mult.append([])
            M.append([])
            M_mult.append([])
            for j in range(N):
                M[i].append([])
                M_mult[i].append([])
#        ret = zeros((N, N))
#        rec_source = zeros(N)
        for i in tqdm(range(N)):
            #%%% Walk through each row (species) %%%
            for dkey, database in self.reactions.items():
                for h123, h123data in database.items():
                    for rkey, r in h123data.items():
                        #print(dkey, h123, rkey)
                        #%%% Sort the species of each reaction into %%%
                        #%%%  the appropriate column %%%
                        # Specify density for external source
                        # Assure that auto-processes are considered
                        # Loop through each reaction defined in the CRM
                        if self.slist[i] in r.educts:
                            M[i][i].append(r.rate)
                            M_mult[i][i].append(-r.e_mult[r.educts.index(self.slist[i])])
                        # Loop through the reaction products
                        for frag in range(len(r.products)):    
                            ''' SOURCE '''
                            # Do nothing if background fragment
                         #   if r.products[frag] not in self.slist: 
                         #       continue
                            # If fragment enters row, add in appropriate column as
                            # defined by reactant
                            if (r.products[frag] in self.slist) and (
                                    self.slist.index(r.products[frag]) == i):
                                # External flag triggered, store to external source
                                common = set(r.educts).intersection(\
                                        set(self.slist))
                                if len(common) == 0:
                                    ''' EXTERNAL SOURCE '''
                                    G[i].append(r.rate)
                                    G_mult[i].append(r.p_mult[frag])
                                # No trigger of external source, store to 
                                # appropriate location in matrix
                                else: 
                                    ''' INTERNAL SOURCE '''
                                    ind = self.slist.index(common.pop())
                                    M[i][ind].append(r.rate)
                                    M_mult[i][ind].append(r.p_mult[frag])
        self.M =  M
        self.M_mult = M_mult
        self.G = G
        self.G_mult = G_mult



    def populate_diag(self, ne, ni=None):
        ''' Function populating and returning the rate matrix and source

        Parameters
        ----------
        ne : float
            electron density in cm**-3
        ni : float, optional (default: None)
            ion density in cm**-3. Default sets ni=ne

        Returns
        -------
        mat : ndarray
            Rate matrix containing the rates at the given 
            temperature and density. A 2D list in case mode=diagnostic
        rec : ndarray
            Vector of recombination sources. List if mode=diagnostic
        '''
        from numpy import zeros,array,sum,transpose
        from math import isinf
        
        if ni is None:
            ni = ne
        N = len(self.species)
        # Setup a 2D diagnostic list for the matrix and a list for 
        # the recombination source
        rec_source = [[] for x in range(N)]
        ret = [[[] for y in range(N)] for x in range(N)]
        for i in range(N):
            #%%% Walk through each row (species) %%%
            for dkey, database in self.reactions.items():
                for h123, h123data in database.items():
                    for rkey, r in h123data.items():
                        #%%% Sort the species of each reaction into %%%
                        #%%%  the appropriate column %%%
                        # TODO: what if three-particle reaction?
                        bg = r.e*ne + r.p*ni # Specify density for reactions
                        # Loop through each reaction defined in the CRM
                        if self.slist[i] in r.educts:
                            ret[i][i].append('-{}*{}_{}{}'.format(
                                    r.e_mult[r.educts.index(self.slist[i])], 
                                    dkey, rkey, bg))
                        # Loop through the reaction products
                        for frag in range(len(r.products)):    
                            ''' SOURCE '''
                            # Do nothing if background fragment
                         #   if r.products[frag] not in self.slist: 
                         #       continue
                            # If fragment enters row, add in appropriate column as
                            # defined by reactant
                            if (r.products[frag] in self.slist) and (
                                    self.slist.index(r.products[frag]) == i):
                                # External flag triggered, store to external source
                                common = set(r.educts).intersection(\
                                        set(self.slist))
                                if len(common) == 0:
                                    ''' EXTERNAL SOURCE '''
                                    ''' Diagnostic matrix '''
                                    rec_source[i].append('+{}*{}_{}{}'.format(
                                            r.p_mult[frag], dkey, rkey, bg))
                                # No trigger of external source, store to 
                                # appropriate location in matrix
                                else: 
                                    ''' INTERNAL SOURCE '''
                                    ''' Diagnostic matrix '''
                                    ret[i][i].append('+{}*{}_{}{}'.format(
                                            r.p_mult[frag], dkey, rkey, bg))
        return ret, rec_source

        
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
        dia, ext = self.populate_diag('*[ne]','*[ni]') 
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

    def getS(self, Te, ne, Ti=None, ni=None, E=0.1,
                                Tm=0,write=False):
        if Ti is None: Ti=Te # Check for Ti, set if necessary
        if ni is None: ni=ne # Check for ni, set if necessary

        mat,ext=self.eval_Sgl(Te=Te,ne=ne,Ti=Ti,ni=ni,Energy=E,Tm=Tm)    


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



      
    def getM(self, Te, ne, Ti=None, ni=None, E=0.1, sparse=False, write=False, 
            ext=None, div=None, extrate=None, mask=True):
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
        ext : N-vector
            external sinks/sources (1/s)
        div : N-vector
            divergence of cell: transport (1/s)
        extrate : N-vector
            external depletion/creation rate (m^3/s)
            
        Returns
        -------
        mat, ext
        mat : ndarray
            Rate matrix containing the rates at the given temperature
            and density in cm**3 s**-1
        ext : ndarray
            Vector of external sinks/sources in s**-1
        '''
        # TODO: make additional sources work with dictionaries not lists

        from scipy.sparse import csc_matrix
        from numpy import zeros, diag

        if Ti is None: 
            Ti = Te # Check for Ti, set if necessary
        if ni is None: 
            ni = ne # Check for ni, set if necessary
        if ext is None:
            ext=zeros(self.N)
        if div is None:
            div=zeros(self.N)
        if extrate is None:
            extrate=zeros(self.N)

        M, rec = self.eval_M(Te=Te, ne=ne, Ti=Ti, ni=ni, E=E)
        #M, rec = self.eval_M(Te, ne, Ti, ni, E)
        if write:   # Write to log if requested
            self.write_matrix(M, rec, 'M', Te, ne, Ti, ne, E)
            if self.verbose: # Print output if running verbose
                with open('logs/M_matrix.log', 'rt') as f:
                    for m in f:
                        print(m.strip())
        if sparse: 
            M = csc_matrix(M) # Use sparse format if requested
        return (M + diag(extrate))*self.evolvemat**mask,\
                (rec + ext + div + self.source)*self.evolvearr**mask
 

    def Sgl(self, Te, ne, Ti=None, ni=None, E=0.1, Tm=0, **kwargs):
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
        **kwargs : keyword arguments
            passed to crm.M
            
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
        mat, ext = self.eval_Sgl(Te=Te, ne=ne, Ti=Ti, ni=ni, E=E, Tm=Tm)    
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
        M, G = self.getM(Te, ne, Ti, ni, E, **kwargs) 
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


    def EI(self, Te, ne, Ti=None, ni=None, E=0.1, Tm=0, write=False):
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
        matE, extE, matI, extI = self.eval_EI(Te=Te, ne=ne, Ti=Ti, ni=ni, Energy=E, Tm=Tm)    
        if write:
            title=['a', 'm']
            for i in range(5):
                self.write_matrix(matI[:,:,i], extI[:,i], 'I'+title[i], Te, ne, Ti,
                        ni, E, form='{:1.2E}')
                self.write_matrix(matE[:,:,i], extE[:,i], 'E'+title[i], Te, ne, Ti,
                        ni, E, form='{:1.2E}')
        return  [[matE[:,:,0], extE[:,0]], [matE[:,:,1], extE[:,1]]], \
                [[matI[:,:,0], extI[:,0]], [matI[:,:,1], extI[:,1]]]
 
    def solve_crm(self, t, Te, ne, Ti=None, ni=None, E=0.1, Tm=0, Srec=True, gl=True,
            n=None, Qres=True, densonly=False,  **kwargs):
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
        Srec : boolean, optional (default: True)
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
        **kwargs : keyword args
            passed to crm.M
        
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
        M, G = self.getM(Te, ne, Ti, ni, E, **kwargs) 
        if gl is True:
            M, G, n0 = self.gl_crm(M, G, Srec, n0)
            sfunc = self.Sgl
        else:
            sfunc = self.getS
        if densonly is False:
            U=sfunc(Te, ne, Ti, ni, E, Tm)
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
                # TODO: Add non-zero ext sources and sinks
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
                # TODO: Add non-zero ext sources and sinks
        else:
            mat = M
            ext = G
            n = n0

        return solve_ivp(lambda x, y: self.ddt(y, x, mat, ext), (0, t),
                n, 'LSODA', dense_output=True)



    def get_transition(self, wavelength, unit='angstrom'):
        from numpy import divide, zeros_like, where
        matE, _, _, _ = self.eval_EI(1,1,1,1,0.1,1)    

        matE = divide(1239.84193, matE, out=zeros_like(matE), where=matE!=0)
        matE = abs(matE-wavelength/(10**(unit=='angstrom')))
        atE = matE[:,:,0]
        molE = matE[:,:,1]
        atl, atu = where(atE == atE.min())
        moll, molu = where(molE == molE.min())
        print('Closest atom transition: {} => {}'.format(self.slist[atu[0]], self.slist[atl[0]]))
        print('Closest molecule transition: {} => {}'.format(self.slist[molu[0]], self.slist[moll[0]]))
        return ((atl[0], atu[0]), (moll[0], molu[0]))
        


    def steady_state(
            self, Te, ne, E=0.1, Ti=None, ni=None, Srec=True, 
            gl=False, plot=False, n0=None, dt=False, store=True,**kwargs):
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
        Srec : boolean, optional (default: True)
            switch for considering 'external' sinks/sources (e.g. 
            re-association)
        gl : boolean, optional (default: False)
            switch determining wheter to evaluate the Np x Np 
            Greenland CRM (True) or solve the full NxN system of ODEs 
            (False)
        plot : boolean, optional (default: False)
            show a plot of the 10s time-evolution
        dt : boolean (default: False)
            Returns the ss solution or the time-trace
        n0 : ndarray, optional (default: None)
            Initial density distribution of the species defined in the 
            input. List must be of the same length as the species list,
            and the order is determined by the input file. Defaults to 
            the initial densities defined in the input if n0 is None.
        **kwargs : keywird arguments
            passed to crm.M

        Returns
        -------
        n_ss : ndarray
            steady-state populations
        '''
        from numpy import log, ones, matmul, polyfit, zeros, exp, sum
        from numpy.linalg import inv
        from scipy.optimize import fsolve
        from scipy.linalg import norm
        from scipy.integrate import solve_ivp
        from matplotlib.pyplot import figure

        if n0 is None:
            n0 = self.n0()
        #Te /= 1.602e-19
        #ne /= 1e6
        # Set volume to be cm**3
        #vol *= 1e6
        # Set initial density in [1/cm**3]
        # Get the rate matrices 
        mat, ext = self.getM(Te, ne, Ti, ni, E, mask=dt,**kwargs) 
        if gl is True:
            # Get matrices and eigenvalues
            mat, ext, n = self.gl_crm(mat, ext, n=n, Srec=Srec) 
        # Set external sources to be [1/cm**3/s]
#        ext[0] = (psorgc + diva + srca + bga)/vol
#        ext[1] = (divm + srcm + bgm)/vol
        # Normalize ionzation rate to UEDGE values 
#        mat[0,0] = (psor/ne)/vol
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
            mat, ext = self.evolve(mat, ext)
            ret = self.fill_static(matmul(inv(mat),-ext))
            if store is True:
                for i in range(self.N):
                    self.species[self.slist[i]]['nss'] = ret[i]
            return ret
            #return fsolve(self.ddt, n0, args=(0, mat, ext))  
        else: 
            ret = solve_ivp(lambda x,y: self.ddt(y, x, mat, ext), (0, 1e10), 
                    n0, 'LSODA', dense_output=True).y[:,-1]
            if store is True:
                for i in range(self.N):
                    self.species[self.slist[i]]['nss'] = ret[i]
            return ret

    def gl_crm(self, te, ne, n0=None, **kwargs):
        ''' Returns the P-space matrices according to Greenland 2001

        Parameters
        ----------
        mat : ndarray
            2D rate matrix
        ext : ndarray
            1D external sink/source vector
        Srec : boolean, optional (default: True)
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
        from numpy import matmul, diag, real, ones
        from numpy.linalg import inv, eig

        # CIRCUMVENT STATIC SPECIES WHICH MAKE MATRICES SINGULAR
        oldarr = self.evolvearr.copy()
        oldmat = self.evolvemat.copy()
        self.evolvearr = ones(self.evolvearr.shape)
        self.evolvemat = ones(self.evolvemat.shape)
        mat, ext = self.getM(te, ne, **kwargs)
        self.evolvearr = oldarr
        self.evolvemat = oldmat 
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
        GPp = real(ext[:self.Np] - matmul(matmul(Delta,
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

    def evaluate_CRM(self, Te, ne, Ti=None, ni=None, E=0.1, printout=True, **kwargs):
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
        **kwargs : keywird arguments
            passed to crm.M

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
        M, T, D = self.gl_crm(*self.getM(Te, ne, Ti, ni, E, **kwargs), 
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

    
    def generate_CRM(self, Te, ne, kappa, Ti=None, ni=None, E=0.1, epsilon=1, **kwargs):
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
        **kwargs : keywird arguments
            passed to crm.M

        Returns
        -------
        None
        '''
        from numpy import array,where,matmul
        from numpy.linalg import inv

        # Get matrices and eigenvalues/-vectors
        M, T, D = self.gl_crm(*self.getM(Te, ne, Ti, ni, E, **kwargs),
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



    def fill_static(self, array):
        from numpy import insert
        for i in range(self.N):
            if self.evolvearr[i] == 0:
                    array = insert(array, i, self.species[self.slist[i]]['n'])
        return array

        
    def intensity(
            self, Te, ne, E=0.1, Ti=None, ni=None, Srec=True, n0=None, 
            units='l', norm=False, write=False, **kwargs):
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
        Srec : boolean, optional (default: True)
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
        from numpy import reshape, zeros, multiply, where, nonzero
        from numpy import transpose, matmul, array, mean, all, diag
        from numpy import divide, zeros_like, where
        from scipy.optimize import minimize
        from scipy.integrate import solve_ivp
        from matplotlib.pyplot import figure


        n_ss = self.steady_state(Te, ne, E=E, Ti=Ti, ni=ni, 
                Srec=Srec, gl=False, plot=False, n0=n0, dt=False, **kwargs)

        x,y = self.EI(Te, ne, Ti, Te, E, write=False)
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
            #arr = transpose(arr[~all(abs(arr) < 1e-8, axis=1)])
            """''' OLD '''
            for nzi in nonzero(arr[:,0]):
                if units == 'v':    
                    arr[nzi,0]=(1e7*arr[nzi,0])/1239.84193
                elif units == 'f':
                    arr[nzi,0]=241.798*arr[nzi,0]
                elif units == 'l':
                    arr[nzi,0]=1239.84193/arr[nzi,0]
                elif units == 'Å':
                    arr[nzi,0]=12398.4193/arr[nzi,0]
            """''' NEW '''
            if units == 'v':    
                arr[:,0]=(1e7*arr[:,0])/1239.84193
            elif units == 'f':
                arr[:,0]=241.798*arr[:,0]
            elif units == 'l':
                arr[:,0] = divide(1239.84193, arr[:,0], out=zeros_like(arr[:,0]), where=arr[:,0]!=0)
            elif units == 'Å':
                arr[:,0] = divide(12398.4193, arr[:,0], out=zeros_like(arr[:,0]), where=arr[:,0]!=0)
            if norm is True:
                nrm = sum(arr[:,1])
            else:
                nrm = 1
            arr[:,1] = arr[:,1]/nrm
            arr = transpose(arr)
            ''' FILTERING START '''
            arr = arr[arr!=0]
            arr = arr.reshape(2,int(len(arr)/2))
            ''' FILTERING END '''
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
            for h123, h123db in database.items():
                for rkey, r in h123db.items():
                    if species in r.educts:
                        if r.type == 'RELAXATION':
                            ret += r.rate(Te, Ti, E=E, ne=ne)/ne
                        else:
                            ret += r.rate(Te, Ti, E=E, ne=ne)
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


    def get_reaction(self, database, h123, name):
        ''' Returns the Reaction object for the requested reaction
        
        Parameters
        ----------
        database : string
            the database where to look for the reaction
        h123 : string
            the reaction type
        name : string
            the reaction name/handle to look for
        '''
         
        return self.reactions[database][h123][name]


    def get_rate(self, database, h123, name, Te, Ti, ne, ni, E=0.1, **kwargs):
        ''' **INCOMPLETE** Returns the rate of the requested reaction '''
        return self.get_reaction(database,h123,name).rate(Te=Te, Ti=Ti, ne=ne, ni=ni, E=E, **kwargs)


    def write_YACORA(self, modelname, path='.', dirname='YACORArun', 
                    Trange=None, write_reac=True, Tlist=None, nlist=None):
        from pathlib import Path
        from numpy import logspace, linspace
        from tqdm import tqdm

        # TODO: spacing, etc

        if Trange is None:
            Trange=logspace(-2, 2.7, 500)
        
        if nlist is None:
            nlist = [1e19, 1e20, 1e21]

        if Tlist is None:
            Tlist = linspace(0.25, 30, 120)

        # Create path for writing    
 #       Path('{}/{}'.format(path, dirname)).mkdir(parents=True, exist_ok=True)
        Path('{}/{}/ratecoefficients'.format(path, dirname)).mkdir(parents=True, exist_ok=True)
 #       Path('{}/{}/modelfiles'.format(path, dirname)).mkdir(parents=True, exist_ok=True)
 #       Path('{}/{}/parameterfiles'.format(path, dirname)).mkdir(parents=True, exist_ok=True)
 #       Path('{}/{}/runfiles'.format(path, dirname)).mkdir(parents=True, exist_ok=True)
 #       Path('{}/{}/output'.format(path, dirname)).mkdir(parents=True, exist_ok=True)


        # Run-file
    
        for n in nlist:
            Path('{}/{}/runfiles/{:.2E}'.format(path, dirname, n)).mkdir(parents=True, exist_ok=True)
            rf = open('{}/{}/runfiles/{:.2E}/run.txt'.format(path, dirname, n),'w')
    
            rf.write('SOLUTION\n')
            rf.write('    NONLINEAR\n\n') 
            
            rf.write('\nMODEL\n')
            rf.write('    ../../modelfiles/{:.2E}/{}.txt\n\n'.format(n, modelname))
            
            rf.write('\nPARAMETERS\n')
            rf.write('    ../../parameterfiles/{:.2E}/param.txt\n\n'.format(n))

            rf.write('\nOUTPUT_FILES\n')
            n_abs = '../../output/{:.2E}/n_abs'.format(n)
            n_bal = '../../output/{:.2E}/n_bal'.format(n)
            pop_coeff = '../../output/{:.2E}/pop_coeff'.format(n)

            Path(n_abs).mkdir(parents=True, exist_ok=True)
            Path(n_bal).mkdir(parents=True, exist_ok=True)
            Path(pop_coeff).mkdir(parents=True, exist_ok=True)
            
            for s in list(self.species)[:-1]:
                rf.write('    {:<24}{:<10}{:<10}{}/{}.txt\n'.format('population_coefficient', s, list(self.species)[0], pop_coeff, s.replace('(','').replace('=','').replace(')','')))
                rf.write('    {:<24}{:<20}{}/{}.txt\n'.format('ABSOLUTE_DENSITY', s, n_abs, s.replace('(','').replace('=','').replace(')','')))
                rf.write('    {:<24}{:<20}{}/{}.txt\n'.format('DENSITY_BALANCE', s, n_bal, s.replace('(','').replace('=','').replace(')','')))
                rf.write('\n')


        # Model file
        for n in nlist:
            Path('{}/{}/modelfiles/{:.2E}'.format(path, dirname, n)).mkdir(parents=True, exist_ok=True)
            mf = open('{}/{}/modelfiles/{:.2E}/{}.txt'.format(path, dirname, n, modelname),'w')
            mf.write('PARTICLES\n')
            for s in self.species:
                mf.write('    {}\n'.format(s))
            mf.write('    H+\n')
            #for bg in self.bg:
            #    mf.write('    {}\n'.format(bg))

            mf.write('\n\nREACTIONS')
            for database, reactions in self.reactions.items():
                for h123, h123db in reactions.items():
                    for handle, rate in h123db.items():
                        mf.write('\n\n')
                        # Write reaction block to file
                        mf.write('    {:<16}{:<16}{:<10}'.format('EDUCTS',rate.educts[1],
                                str(int(rate.e_mult[1]))))
                        if rate.educts[0] == 'e':
                            filefac = 'NE'
                        else:
                            filefac = 'NI1'
                        mf.write('\n')

                        mf.write('    {:<16}'.format('PRODUCTS'))
                        for fi in range(len(rate.products)):
                            species = rate.products[fi].replace('p', 'H+')
                            if rate.products[fi] == rate.educts[0]:
                                if rate.products[fi] == 'e':
                                    pass
                                elif rate.p_mult[fi] > rate.e_mult[0]:
                                    mf.write('{:<16}{:<10}'.format(species, str(int(rate.p_mult[fi] - rate.e_mult[0]))))
                                else: 
                                    pass
                            else:
                                mf.write('{:<16}{:<10}'.format(species, str(int(rate.p_mult[fi]))))
                        mf.write('\n')                

                        if rate.type == 'H.4':
                            ratepath =  '../../ratecoefficients/{:.2E}/{}-{}.txt'.format(n,database, handle)
                            commentstr = '{} rate coefficient for reaction {} at ne={:.2E}'.format(database, handle, n)
                        else:
                            ratepath = '../../ratecoefficients/{}-{}.txt'.format(database, handle)
                            commentstr = '{} rate coefficient for reaction {}'.format(database, handle)
                        mf.write('    {:<16}{:<16}{:<56}    {}\n'.format('PROBABILITY','FILEFACTORS', ratepath, filefac))
                        mf.write('    {:<16}{}'.format('COMMENT', commentstr))

                        if write_reac is True:
                            dbpath='{}/{}/ratecoefficients'.format(path, dirname)
                            if rate.type == 'H.4': # Density-dependent rate
                                ndbpath = '{}/{:.2E}'.format(dbpath,n)
                                Path(ndbpath).mkdir(parents=True, exist_ok=True)
                                f = open('{}/{}-{}.txt'.format(ndbpath, database, 
                                    handle),'w')
                                f.write('RATE_COEFFICIENT\n')
                                f.write('VALUES\n')
                                f.write((8*' ').join(['Te','X\n']))
                                for T in Trange:
                                    f.write('{:.5E}  {:.5E}\n'.format(T, 
                                        1e-6*rate.rate(Te=T, Ti=T, ne=n*1e-6, ni=n*1e-6, E=0.1)))
                            elif rate.type == 'H.3': # Energy-dependent rate
                                f = open('{}/{}-{}.txt'.format(dbpath, database, 
                                    handle),'w')
                                f.write('RATE_COEFFICIENT\n')
                                f.write('VALUES\n')
                                f.write((8*' ').join(['Te','X\n']))
                                for T in Trange:
                                    f.write('{:.5E}  {:.5E}\n'.format(T, 
                                        rate.rate(T, T, E=0.1)))
                            else:
                                f = open('{}/{}-{}.txt'.format(dbpath, database, 
                                    handle),'w')
                                f.write('RATE_COEFFICIENT\n')
                                f.write('VALUES\n')
                                f.write((8*' ').join(['Te','X\n']))
                                for T in Trange:
                                    f.write('{:.5E}  {:.5E}\n'.format(T, 
                                        1e-6*rate.rate(Te=T, Ti=T, ne=n*1e-6, ni=n*1e-6, E=0.1)))

        # Runfile
        e2k = 1.160451812e4
        for n in nlist:
            Path('{}/{}/parameterfiles/{:.2E}'.format(path, dirname, n)).mkdir(parents=True, exist_ok=True)
            pf = open('{}/{}/parameterfiles/{:.2E}/param.txt'.format(path, dirname, n),'w')
            pf.write('PROFILE\n')
            pf.write('    TE\n    TI1\n\n')

            pf.write('\nTEMPERATURES\n')
            pf.write('TE\n    VALUES\n')
            pf.write('    {}\n'.format(int(len(Tlist))))
            for T in Tlist:
                pf.write('    {:.2E}\n'.format(T))
            pf.write('\nTI1\n    VALUES\n')
            pf.write('    {}\n'.format(int(len(Tlist))))
            for T in Tlist:
                pf.write('    {:.2E}\n'.format(T*e2k))
            pf.write('\nTN1\n    FIXED\n    {:.2E}\n\n'.format((2/5)*0.1*e2k)) # E=0.1eV

            pf.write('\nDENSITIES\n')
            pf.write('NE\n    FIXED\n    {:.2E}\n\n'.format(n))
            pf.write('NI1\n   FIXED\n    {:.2E}\n\n'.format(n))
            pf.write('NN1\n   FIXED\n    1e18\n\n')
            
            pf.write('\nSTARTING_DENSITIES\n')
            pf.write('    {:<14}{:<10}FIXED\n'.format(list(self.species)[0], 'NN1'))
            for i in range(1,15):                       
                pf.write('    {:<14}{:<10}VARIABLE\n'.format(list(self.species)[i], '0'))
            pf.write('    {:<14}{:<10}FIXED\n'.format(list(self.species)[15], '0'))
            pf.write('    {:<14}{:<10}FIXED\n'.format('H+', 'NI1'))

            pf.write('\n\nEEPF\n    MAXWELL')
            
                



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

    def evolve(self, mat, rec):
        # TODO: fix steady-state when two or more static species!
        from numpy import delete

        for i in range(len(self.evolvearr)-1,-1,-1):
            if self.evolvearr[i] == 0:
                mat = delete(mat, i, axis=0)
                rec = delete(rec, i, axis=0)
                rec += mat[:,i]*self.n0()[i]
                mat = delete(mat, i, axis=1)

        return mat, rec












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
                for rea in range(len(r.educts)):
                    # Find the column into which the products goes: if background mark external
                    try:
                        j=self.slist.index(r.educts[rea])  # Get the product species index of the correct column
                    except:
                        continue
            
                    multiplier=r.e_mult[rea] # Get the reaction multiplier for the density

                    if self.slist[i]==r.educts[rea]:   # If the species (row index) is a reactant, r is a depletion process 
                        ''' DEPLETION '''

                        if mode=='diagnostic':
                            ''' Diagnostic matrix '''
                            ret[i][i].append('-'+str(multiplier)+'*'+r.database+'_'+r.name+bg)   # Print the rate to the correct element

                        elif mode=='R':
                            ''' Rate coefficient matrix '''
                            ret[i,i]-=multiplier*r.rate(Te,Ti,E=E,ne=ne) # Calculate the Rate coefficient and store appropriately

                        elif mode=='M':
                            ''' Rate matrix '''
                            ret[i,i]-=multiplier*r.rate(Te,Ti,E=E,ne=ne)*bg # Calculate the rate and store appropriately

                for frag in range(len(r.products)):    # Loop through the reaction products
                    ''' SOURCE '''
                    multiplier=r.p_mult[frag]   # Fragment multiplier

                    # Do nothing if background fragment
                    if r.products[frag] not in self.slist: 
                        continue
                    
                    # If fragment enters row, add in appropriate column as defined by reactant
                    elif self.slist.index(r.products[frag])==i:
                            if j is None: # External flag triggered, store to external source
                                ''' EXTERNAL SOURCE '''

                                if mode=='diagnostic':
                                    ''' Diagnostic matrix '''
                                    ext_source[i].append('+'+str(multiplier)+'*'+r.database+'_'+r.name+bg)

                                elif mode=='R':
                                    ''' Rate coefficient matrix '''
                                    ext_source[i]+=multiplier*r.rate(Te,Ti,E=E,ne=ne)

                                elif mode=='M':
                                    ''' Rate matrix '''
                                    ext_source[i]+=multiplier*r.rate(Te,Ti,E=E,ne=ne)*bgm

                                elif mode=='Sgl':
                                    ''' Energy source matrix in Greenland form '''
                                    ext_source[i,:]+=r.rate(Te,Ti,E=E,ne=ne)*bgm*Sgl[:,0]+Sgl[:,1]
            
                                elif mode=='I':
                                    ''' Intensity matrix '''
                                    ext_source[i,:]+=r.rate(Te,Ti,E=E,ne=ne)*bgm*Sgl[:,0]

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
                                    ret[i,j]+=multiplier*r.rate(Te,Ti,E=E,ne=ne)

                                elif mode=='M':
                                    ''' Rate matrix '''
                                    ret[i,j]+=multiplier*r.rate(Te,Ti,E=E,ne=ne)*bg

                                elif mode=='Sgl': 
                                    ''' Energy source matrix in Greenland form '''
                                    ret[i,j,:]+=r.rate(Te,Ti,E=E,ne=ne)*bg*Sgl[:,0]+Sgl[:,1]

                                elif mode=='I':
                                    ''' Correlation matrix '''    
                                    ret[i,j,:]+=r.rate(Te,Ti,E=E,ne=ne)*bg*(abs(Sgl[:,0])>0)

                                elif mode=='E':
                                    ''' Correlation matrix '''    
                                    ret[i,j,:]+=Sgl[:,0]
                                    
        return ret,ext_source



    def full_nt(self,Te,ne,t,gl,Ti=None,ni=None,E=0.1,n=None,Srec=True, **kwargs):
        ''' **OBSOLETE** '''
        from scipy.integrate import solve_ivp 

        # Use input n0 as default unless explict n0 requested
        if n is None: n=self.n0() 

        if gl is True:
            # Set up Greenland model 
            mat,ext,nP0p=self.gl_crm(*self.getM(Te,ne,Ti,ni,E,**kwargs),Srec,n)
            n=nP0p
        else:
            mat,ext=self.getM(Te,ne,Ti,ni,E,**kwargs) # Get the full rate matrix

        # Solve and return
        return solve_ivp(lambda x,y: self.solve_crm(x,y,mat,ext),
                                            (0,t),n,method='LSODA') 



        
    def dEdt_old(self,t,Te,ne,Ti=None,ni=None,E=0.1,Tm=0,
                Srec=True,write=False,gl=True,n=None,Qres=True):
        ''' **OBSOLETE** '''
        from numpy import block,zeros,matmul,reshape,sum
        from  numpy.linalg import inv
        from scipy.integrate import solve_ivp


        N=len(self.species)
        Np=self.Np
        Nq=N-Np  


        

        M,G=self.getM(Te,ne,Ti,ni,E,write=write) # Get the full rate matrix

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

            Meff,GPp,nP0p=self.gl_crm(M,G,Srec,n0)
        
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

        


    def populate_new_compound(self, mode, Te, ne, Ti=None, ni=None, E=0.1, Tm=0):
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
        rec : ndarray
            Vector of recombination sources. List if mode=diagnostic
        '''
        from numpy import zeros,array,sum,transpose
        from math import isinf
        
        N = len(self.species)
        if mode == 'diagnostic':
            # Setup a 2D diagnostic list for the matrix and a list for 
            # the recombination source
            rec_source = [[] for x in range(N)]
            ret = [[[] for y in range(N)] for x in range(N)]
        elif mode in ['E', 'I', 'Sgl']:
            ret = zeros((N, N, 2 + 3*(mode == 'Sgl')))
            rec_source = zeros((N, 2+3*(mode == 'Sgl')))
        else:
            # Setup a matrix and vector
            ret = zeros((N, N))
            rec_source = zeros(N)
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
                        rate = r.rate(Te, Ti, E=E, ne=ne)
                        # Fix for writing AMJUEL 2.2.9 for high temperatures 
                        # and density - required for UEDGE fits but
                        # outside the polynomial fit range
                        if isinf(rate):
                            rate = 1e20
                    if mode == 'R':
                        bgm = rate
                        bg = rate
                    elif mode in ['M']:#,'Sgl']:
                        # Specify density for external source
                        bgm = bgm*rate
                        # Assure that auto-processes are considered
                        bg = max(bg, 1)*rate
                    elif mode == 'Sgl':
                        # Specify density for external source
                        bgm = bgm*rate*Sgl[:,0] 
                        # Assure that auto-processes are considered
                        bg = max(bg, 1)*rate*Sgl[:,0]
                        ext = Sgl[:,1]
                    elif mode == 'I':
                        print(rate)
                        # Specify density for external source
                        bgm = bgm*rate*Sgl[-2:,0] 
                        # Assure that auto-processes are considered
                        bg = max(bg, 1)*rate*(abs(Sgl[-2:,0]) > 0)    
                    elif mode == 'E':
                        bgm = Sgl[-2:,0] # Specify density for external source
                        # Assure that auto-processes are considered
                        bg = Sgl[-2:,0]  
                    j = None # Set flag to identify external sources
                    # Loop through each reaction defined in the CRM
                    for rea in range(len(r.educts)):
                        # Find the column into which the products goes: 
                        # if background mark external
                        try:
                            # Get the product species index of the 
                            # correct column
                            j = self.slist.index(r.educts[rea])  
                        except:
                            continue
                        # If the species (row index) is a reactant, r is a 
                        # depletion process 
                        if self.slist[i] == r.educts[rea]:   
                            ''' DEPLETION '''
                            if mode == 'diagnostic':
                                ''' Diagnostic matrix '''
                                # Print the rate to the correct element
                                ret[i][i].append('-{}*{}_{}{}'.format(
                                        r.e_mult[rea], dkey, rkey, bg))
                            elif mode in ['R','M']:
                                ''' Rate coefficient matrix '''
                                # Calculate the Rate coefficient and store 
                                # appropriately
                                ret[i,i]-=r.e_mult[rea]*bg 
                    # Loop through the reaction products
                    for frag in range(len(r.products)):    
                        ''' SOURCE '''
                        # Do nothing if background fragment
                     #   if r.products[frag] not in self.slist: 
                     #       continue
                        # If fragment enters row, add in appropriate column as
                        # defined by reactant
                        if (r.products[frag] in self.slist) and (
                                self.slist.index(r.products[frag]) == i):
                            multiplier = r.p_mult[frag]**(mode not in 
                                    ['Sgl', 'I', 'E'])   # Fragment multiplier
                            # External flag triggered, store to external source
                            if j is None: 
                                ''' EXTERNAL SOURCE '''
                                if mode == 'diagnostic':
                                    ''' Diagnostic matrix '''
                                    rec_source[i].append('+{}*{}_{}{}'.format(
                                            multiplier, dkey, rkey, bg))
                                else:
                                    try:
                                        rec_source[i] += multiplier*bgm + ext
                                    except:
                                        rec_source[i,:] += multiplier*bgm + ext
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
        return ret, rec_source

      
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
        R, ext = self.eval_R(Te=Te, ne=0, Ti=Ti, ni=0, E=E)
        if write: # Write to log if requested
            self.write_matrix(R, ext, 'R', Te, 0, Ti, 0, E)
            if self.verbose: # Print rate matrix to stdout if running verbose
                with open('logs/R_matrix.log', 'rt') as f:
                    for m in f:
                        print(m.strip())
        if sparse is True: 
            R = csc_matrix(R)  # Use sparse format if requested
        return R, ext


    # TODO: make create_Sgl function and eval_Sgl function
    def populate_Sgl(self, Te, ne, Ti=None, ni=None, E=0.1, Tm=None):
        ''' Function populating and returning the rate matrix and source

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
        Tm : float, optional (default: None)
            local molecular temperature for calculating energy losses 
            associated with molecules. Defaults to 2E/3 if Tm=None.

        Returns
        -------
        mat : ndarray
            Rate matrix containing the rates at the given 
            temperature and density. A 2D list in case mode=diagnostic
        rec : ndarray
            Vector of recombination sources. List if mode=diagnostic
        '''
        from numpy import zeros,array,sum,transpose
        from math import isinf
        
        if Tm is None:
            Tm = (2/3)*E
        if ni is None:
            ni = ne
        if Ti is None:
            Ti = Te
        N = len(self.species)
        ret = zeros((N, N, 5))
        rec_source = zeros((N, 5))
        count=0
        for i in range(N):
            #%%% Walk through each row (species) %%%
            for dkey, database in self.reactions.items():
                for h123, h123data in database.items():
                    for rkey, r in h123data.items():
                        #%%% Sort the species of each reaction into %%%
                        #%%%  the appropriate column %%%
                        Sgl = r.getS(Te, Ti, Tm, E)
                        # TODO: what if three-particle reaction?
                        rate = r.rate(Te=Te, Ti=Ti, E=E, ne=ne, ni=ni)
                        # Fix for writing AMJUEL 2.2.9 for high temperatures 
                        # and density - required for UEDGE fits but
                        # outside the polynomial fit range
                        if isinf(rate):
                            rate = 1e20
                        # Loop through the reaction products
                        for frag in range(len(r.products)):    
                            ''' SOURCE '''
                            # Do nothing if background fragment
                         #   if r.products[frag] not in self.slist: 
                         #       continue
                            # If fragment enters row, add in appropriate column as
                            # defined by reactant
                            if (r.products[frag] in self.slist) and (
                                    self.slist.index(r.products[frag]) == i):
                                # External flag triggered, store to external source
                                common = set(r.educts).intersection(\
                                        set(self.slist))
                                if len(common) == 0:
                                    ''' EXTERNAL SOURCE '''
                                    rec_source[i,:] += ( \
                                        (r.e*ne) * (r.p*ni))*rate*Sgl[:,0] + Sgl[:,1]
                                # No trigger of external source, store to 
                                # appropriate location in matrix
                                else: 
                                    ''' INTERNAL SOURCE '''
                                    count += 1
                                    ind = self.slist.index(common.pop())
                                    ret[i,ind,:] += \
                                        max((r.e*ne + r.p*ni), 1)*rate*Sgl[:,0] \
                                        + Sgl[:,1]
        return ret, rec_source


    def populate_R(self, Te, ne, Ti=None, ni=None, E=0.1):
        ''' Function populating and returning the rate matrix and source

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

        Returns
        -------
        mat : ndarray
            Rate matrix containing the rates at the given 
            temperature and density. A 2D list in case mode=diagnostic
        rec : ndarray
            Vector of recombination sources. List if mode=diagnostic
        '''
        from numpy import zeros,array,sum,transpose
        from math import isinf
        
        if ni is None:
            ni = ne
        if Ti is None:
            Ti = Te
        N = len(self.species)
        # Setup a matrix and vector
        ret = zeros((N, N))
        rec_source = zeros(N)
        for i in range(N):
            #%%% Walk through each row (species) %%%
            for dkey, database in self.reactions.items():
                for h123, h123data in database.items():
                    for rkey, r in h123data.items():
                        #print(dkey, h123, rkey)
                        #%%% Sort the species of each reaction into %%%
                        #%%%  the appropriate column %%%
                        # TODO: what if three-particle reaction?
                        rate = r.rate(Te, Ti, E=E, ne=ne, ni=ni)
                        # Fix for writing AMJUEL 2.2.9 for high temperatures 
                        # and density - required for UEDGE fits but
                        # outside the polynomial fit range
                        if isinf(rate):
                            rate = 1e20
                        # Specify density for external source
                        # Assure that auto-processes are considered
                        # Loop through each reaction defined in the CRM
                        if self.slist[i] in r.educts:
                            ret[i,i]-=r.e_mult[r.educts.index(self.slist[i])]*rate
                        # Loop through the reaction products
                        for frag in range(len(r.products)):    
                            ''' SOURCE '''
                            # Do nothing if background fragment
                         #   if r.products[frag] not in self.slist: 
                         #       continue
                            # If fragment enters row, add in appropriate column as
                            # defined by reactant
                            if (r.products[frag] in self.slist) and (
                                    self.slist.index(r.products[frag]) == i):
                                # External flag triggered, store to external source
                                common = set(r.educts).intersection(\
                                        set(self.slist))
                                if len(common) == 0:
                                    ''' EXTERNAL SOURCE '''
                                    rec_source[i] += r.p_mult[frag]*rate
                                # No trigger of external source, store to 
                                # appropriate location in matrix
                                else: 
                                    ''' INTERNAL SOURCE '''
                                    ret[i,self.slist.index(common.pop())] += \
                                        r.p_mult[frag]*rate

        return ret, rec_source





    def populate_M(self, Te, ne, Ti=None, ni=None, E=0.1):
        ''' Function populating and returning the rate matrix and source

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

        Returns
        -------
        mat : ndarray
            Rate matrix containing the rates at the given 
            temperature and density. A 2D list in case mode=diagnostic
        rec : ndarray
            Vector of recombination sources. List if mode=diagnostic
        '''
        from numpy import zeros,array,sum,transpose
        from math import isinf
        
        if ni is None:
            ni = ne
        if Ti is None:
            Ti = Te
        N = len(self.species)
        # Setup a matrix and vector
        ret = zeros((N, N))
        rec_source = zeros(N)
        for i in range(N):
            #%%% Walk through each row (species) %%%
            for dkey, database in self.reactions.items():
                for h123, h123data in database.items():
                    for rkey, r in h123data.items():
                        #print(dkey, h123, rkey)
                        #%%% Sort the species of each reaction into %%%
                        #%%%  the appropriate column %%%
                        # TODO: what if three-particle reaction?
                        rate = r.rate(Te, Ti, E=E, ne=ne, ni=ni)
                        # Fix for writing AMJUEL 2.2.9 for high temperatures 
                        # and density - required for UEDGE fits but
                        # outside the polynomial fit range
                        if isinf(rate):
                            rate = 1e20
                        # Specify density for external source
                        bg = ((r.e*ne) * (r.p*ni))
                        # Assure that auto-processes are considered
                        # Loop through each reaction defined in the CRM
                        if self.slist[i] in r.educts:
                            ret[i,i]-=r.e_mult[r.educts.index(self.slist[i])]*max(bg, 1)*rate
                        # Loop through the reaction products
                        for frag in range(len(r.products)):    
                            ''' SOURCE '''
                            # Do nothing if background fragment
                         #   if r.products[frag] not in self.slist: 
                         #       continue
                            # If fragment enters row, add in appropriate column as
                            # defined by reactant
                            if (r.products[frag] in self.slist) and (
                                    self.slist.index(r.products[frag]) == i):
                                # External flag triggered, store to external source
                                common = set(r.educts).intersection(\
                                        set(self.slist))
                                if len(common) == 0:
                                    ''' EXTERNAL SOURCE '''
                                    rec_source[i] += r.p_mult[frag]*bg*rate
                                # No trigger of external source, store to 
                                # appropriate location in matrix
                                else: 
                                    ''' INTERNAL SOURCE '''
                                    ret[i,self.slist.index(common.pop())] += \
                                        r.p_mult[frag]*max(bg, 1)*rate

        return ret, rec_source


    # TODO: make create_EI function and eval_EI function
    def populate_EI(self, Te, ne, Ti=None, ni=None, E=0.1, Tm=None):
        ''' Function populating and returning the rate matrix and source

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
        Tm : float, optional (default: 0)
            local molecular temperature for calculating energy losses 
            associated with molecules. Defaults to E if Tm=0.

        Returns
        -------
        matE : ndarray
            Rate matrix containing the rates at the given 
            temperature and density. A 2D list in case mode=diagnostic
        recE : ndarray
            Vector of recombination sources. List if mode=diagnostic
        matI : ndarray
            Rate matrix containing the rates at the given 
            temperature and density. A 2D list in case mode=diagnostic
        recI : ndarray
            Vector of recombination sources. List if mode=diagnostic
        '''
        from numpy import zeros,array,sum,transpose
        from math import isinf
        
        if Tm is None:
            Tm = E
        if ni is None:
            ni = ne
        if Ti is None:
            Ti = Te
        N = len(self.species)
        retE = zeros((N, N, 2))
        rec_sourceE = zeros((N, 2))
        retI = zeros((N, N, 2))
        rec_sourceI = zeros((N, 2))
        for i in range(N):
            #%%% Walk through each row (species) %%%
            for dkey, database in self.reactions.items():
                for h123, h123data in database.items():
                    for rkey, r in h123data.items():
                        #%%% Sort the species of each reaction into %%%
                        #%%%  the appropriate column %%%
                        Sgl = r.getS(Te, Ti, Tm, E)
                        # TODO: what if three-particle reaction?
                        # Fix for writing AMJUEL 2.2.9 for high temperatures 
                        # and density - required for UEDGE fits but
                        # outside the polynomial fit range
                        rate = r.rate(Te, Ti, E=E, ne=ne, ni=ni)
                        if isinf(rate):
                            rate = 1e20
                        for frag in range(len(r.products)):    
                            ''' SOURCE '''
                            # Do nothing if background fragment
                         #   if r.products[frag] not in self.slist: 
                         #       continue
                            # If fragment enters row, add in appropriate column as
                            # defined by reactant
                            if (r.products[frag] in self.slist) and (
                                    self.slist.index(r.products[frag]) == i):
                                # External flag triggered, store to external source
                                common = set(r.educts).intersection(\
                                        set(self.slist))
                                if len(common) == 0:
                                    ''' EXTERNAL SOURCE '''
                                    rec_sourceI[i,:] += ((r.e*ne) * \
                                        (r.p*ni))*rate*Sgl[-2:,0] 
                                    rec_sourceE[i,:] += Sgl[-2:,0]
                                # No trigger of external source, store to 
                                # appropriate location in matrix
                                else: 
                                    ''' INTERNAL SOURCE '''
                                    j = self.slist.index(common.pop())
                                    # TODO: Other check than positive energy change?
                                    # How should "spontaneous excitation" be considered?
                                    retI[i,j,:] += max(r.e*ne + r.p*ni, 1)*rate* \
                                            (Sgl[-2:,0] > 0)   
                                    retE[i,j,:] += Sgl[-2:,0]*(Sgl[-2:,0] >0)
        return retE, rec_sourceE, retI, rec_sourceI



    def old_steady_state(
            self, Te, ne, na, nm, diva, divm, srca, srcm, bga, bgm, psor,
            psorgc, vol, E=0.1, Ti=None, ni=None, Srec=True, gl=False,
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
        Srec : boolean, optional (default: True)
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
        mat, ext = self.getM(Te, ne, Ti, ni, E, write=True) 
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


