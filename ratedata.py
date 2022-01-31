# Class for reading and storing rate data from different locations
# Created based on read_EIRENE.py by holm10 on Jan 27 2020
# Changelog:
#   200127 - Rewrote read_EIRENE into a class


class RateData:
    ''' Class containing data stored in the default databases

    This class sets up an atomic and molecular reaction rate database
    from a set of default rate databases. 

    ...

    Attributes
    ----------
    reactions : dict
        dictionary containing a sub-dictionary for each default database
        'AMJUEL' : dict
            Handles for each reaction, containing the rate coefficient 
            from the AMJUEL database
        'HYDHEL' : dict
            Handles for each reaction, containing the rate coefficient 
            from the HYDHEL database
        'H2VIBR' : dict
            Handles for each reaction, containing the rate coefficient 
            from the H2VIBR database
        'ADAS' : dict
            Handles for each reaction, containing the rate coefficient 
            from the ADAS database
        'UE' : dict
            Handles for each reaction, containing the rate coefficient 
            from the UE database
    '''


    def __init__(
            self, rates = {'AMJUEL':None, 'HYDHEL':None, 'H2VIBR':None, 
            'ADAS':None, 'UE':None}, path='.'):
        ''' 
        Parameters
        ----------
        amjuel : string, optional (default: 'amjuel.tex')
            path to AMJUEL LaTeX file, relative to path 
        hydhel : string, optional (default: 'hydhel.tex')
            path to HYDHEL LaTeX file, relative to Crm.path 
        h2vibr : string, optional (default: 'h2vibr.tex')
            path to H2VIBR LaTeX file, relative to Crm.path 
        ADAS : string, optional (default: 'ich0-1.dat')
            path to ADAS LaTeX file, relative to Crm.path 
        UE : string, optional (default: 'ehr1.dat')
            path to UEDGE DEGAS2-style LaTeX file, relative to Crm.path
        path : string, optional (default: '.')
            path, to which all other paths are relative
        '''
        # Create a loop that reads the EIRENE tex files, compatible with
        # Jan 2020 versions
        self.reactions = {'FCF':{},'AIK':{}}
        for database, subpath  in rates.items():
            if (database.upper() == 'AMJUEL') and (subpath is not None):
                self.reactions[database.upper()] = {}
                self.read_EIRENE(subpath, self.reactions[database], 
                        [500, ['b0','0','a0','h0','p0','k0'], ['a0','h0',
                        'p0','k0'], 45],  path=path)
            elif (database.upper() == 'HYDHEL') and (subpath is not None):
                self.reactions[database.upper()] = {}
                self.read_EIRENE(subpath, self.reactions[database], [150, 
                        ['b0','0','a0','h0'], ['a0','h0'], 80], path=path)
            elif (database.upper() == 'H2VIBR') and (subpath is not None):
                self.reactions[database.upper()] = {}
                self.read_EIRENE(subpath, self.reactions[database], [0, 
                        ['b0','0','a0','h0'], ['a0','h0'], 20], path=path)
            elif database[:3].upper() == 'FCF':
                self.reactions['FCF'][database.split('_')[1].strip()] = \
                        self.read_factors(subpath)
            elif database[:3].upper() == 'AIK':
                self.reactions['AIK'][database.split('_')[1].strip()] = \
                        self.read_factors(subpath)
                

        '''
        if ('AMJUEL' in list(rates)) and (rates['AMJUEL'] is not None):
            self.reactions['AMJUEL'] = { 'settings' : [500, 
                            ['b0','0','a0','h0','p0','k0'],
                            ['a0','h0','p0','k0'], 45], 
                            'path':rates['AMJUEL']}
        if ('HYDHEL' in list(rates)) and (rates['HYDHEL'] is not None):
            self.reactions['HYDHEL'] = {'settings' : [150, 
                            ['b0','0','a0','h0'], ['a0','h0'], 80], 
                            'path' : rates['HYDHEL'] }
        if ('H2VIBR' in list(rates)) and (rates['H2VIBR'] is not None):
            self.reactions['H2VIBR'] = {'settings' : [0, ['b0','0','a0','h0'],
                            ['a0','h0'], 20], 'path' : rates['H2VIBR']}

        # For each data point, add the reactions to the appropriate dictionary
        for rate in list(rates):
            try:
                self.read_EIRENE(self.reactions[rate]['path'],
                        self.reactions[rate], 
                        self.reactions[rate]['settings'], path=path)
            except:
                print(  'Database {} not found in {}/{}. Omitting'
                        ''.format(rate, path, 
                        self.reactions[rate]['path']))
        if ('ADAS' in list(rates)) and (rates['ADAS'] is not None):
            self.reactions['ADAS'] = {}
            try:
                self.read_ADAS(rates['ADAS'], self.reactions['ADAS'], path=path)
            except:
                print(  'Database ADAS not found in {}/{}. Omitting'
                        ''.format(path, ADAS))
        if ('UE' in list(rates)) and (rates['UE'] is not None): 
            self.reactions['UE'] = {}
            try:
                self.read_UE(rates['UE'],self.reactions['UE'],path=path)
            except:
                print(  'Database UE not found in {}/{}. Omitting'
                        ''.format(path, UE)) '''


    def read_factors(self, fname):
        from numpy import array, transpose
        data = []
        with open(fname, 'r') as f:
            i = 0
            for line in f:
                if i>0:
                    data.append([float(x) for x in line.split()])
                i += 1
        data = array(data)

        return transpose(array(data))[1:]
        

    def read_EIRENE(self, fname, reactions, settings, path='.'):
        ''' Reads reaction coefficents from EIRENE databases

        Created to read EIRENE AMJUEL, HYDHEL, and H2VIBR data from
        the LaTeX files available from eirene.de

        fname : string
            file name of file to be read
        reactions : dict
            location where the reaction handles are saved to
        settings : list
            list of settings specifying the behaviour of the read 
            routine, with the four following elements
            settings[0] : int
                lines to omit from start of file
            settings[1] : list of strings
                non-coefficient rate entries: omitted
            settings[2] : list of strings 
                data rate entries to be read (?)
            settings[3] : int
                lines to omit from end of file
        path : string, optional (default: '.')
            path relative to which fname is searched for

        Returns
        -------
        None
        '''
        from numpy import array
        # Script parsing Hydhel.tex
        book = []

        data = False
        # Open and loop through file line-by-line
        with open('{}/{}'.format(path, fname), 'r') as f:
            for line in f:
                # Find the beginning of data
                if '##BEGIN DATA HERE##' in line:
                    data = True
                elif 'section{Appendix}' in line:
                    data = False
                elif data is False:
                    pass
                else:
                    if len(line.split()) !=0:
                        book.append(line.split())
        h123 = None

        oneparfit = ['a', 'b', 'e', 'h', 'k', 'p']
        twoparfit = ['H.3', 'H.4', 'H.7', 'H.9', 'H.10', 'H.12']
        reaction = False
        coeffs = False
        for l in book:
            # Reaction is being read
            if reaction is True:
                # Check whether we're at the end of coefficient declatation
                if 'end{verbatim' in l[0]:
                    reaction = False
                    coeffs = False
                    reactions[h123][rname.upper()] = array(cf)
                # Double polynomial fit
                if h123 in twoparfit:
                    if 'begin{verbatim' in l[0]:
                        coeffs = True
                    elif coeffs is True and l[0] in [str(x) for x in range(9)]:
                        cf[int(l[0])] += [float(x.replace('D','E')) for x in l[1:]]
                # Single polynomial fit
                else:
                    if 'begin{verbatim' in l[0]:
                        coeffs = True
                    elif coeffs is True:
                        if l[0][0] in oneparfit:
                            cf += [float(x.strip(',').replace('D','E')) for x in l[1::2]]
            # Begin reading reaction
            elif l[0] == 'Reaction':
                reaction = True # Reaction coeffs follow
                rname = l[1] # Store the reaction name
                if h123 in twoparfit:
                    cf = [[] for x in range(9)]
                else:
                    cf = []
            # Move to other fit type
            elif r'\section' in l[0]:
                h123 = l[0].split('{')[1].strip(':')#.replace('.','')
                reactions[h123] = dict()

    def read_ADAS(self, fname, reactions, path='.'):
        ''' Reads coefficents from ADAS databases

        Created to read ADAS data from OPEN-ADAS data files

        fname : string
            file name of file to be read
        reactions : dict
            location where the reaction handles are saved to
        path : string, optional (default: '.')
            path relative to which fname is searched for

        Returns
        -------
        None
        '''
        # Constants turning the values into 
        cm1 = 1/8065.6
        kB = 8.621738e-5
        # The ADAS reading routine follows the ADAS documentation: only 
        # parameters relevant to the CRM are extracted
        with open('{}/{}'.format(path, fname)) as f:  # Open the file
            l = f.readline() # Discard data
            l = f.readline().split() # Read first data line
            # Loop through energy level data
            while l[0] != '-1':
                reactions['E' + l[0].upper()] = float(l[4])*cm1 # Store as eV
                                                        # in reactions['E']
                l = f.readline().split()
            # Read and store temperature point data
            l = f.readline().split()
            reactions['T'] = [float(x.replace('+', 'e+').replace(
                    '-', 'e-'))*kB for x in l[2:]]
            # Read first rate line
            l = f.readline()
            # Read all unspecified data
            while l[0] == ' ':
                if l.strip() == '-1': break
                [ul, ll] = l[:8].split()
                l = l[8:]
                l = [l[i:i+8].strip() for i in range(0, len(l), 8)][:-1]
                # Store excitation and relaxation data
                reactions[ul.upper() + '-' + ll.upper()] = float(l[0][:-3] \
                        + 'e' + l[0][-3:])
                reactions[ll.upper() + '-' + ul.upper()] = [float(x[:-3] \
                        + 'e' + x[-3:]) for x in l[1:len(reactions['T']) + 1]]
                l = f.readline()


    def read_UE(self, fname, reactions, path='.',
                datalist=['IONIZ', 'REC', 'IONIZRAD', 'RECRAD']):
        ''' Reads coefficents from ADAS databases

        Created to read ADAS data from OPEN-ADAS data files

        fname : string
            file name of file to be read
        reactions : dict
            location where the reaction handles are saved to
        path : string, optional (default: '.')
            path relative to which fname is searched for
        datalist : list of strings, optional 
                    (default: ['IONIZ', 'REC', 'IONIZRAD', 'RECRAD'])
            List of names to give to blocks in rate data file. All 
            blocks are read in consecutive order and stores with the 
            corresponding datalist entry as their database.
        Returns
        -------
        None
        '''
        from numpy import array, transpose

        with open('{}/{}'.format(path, fname)) as f:  # Open the file
            # Read the ionizationa nd recombination data (1st and 2nd blocks
            # in ehr1.dat)
            for k in datalist:
                k = k.upper()
                reactions[k] = []
                l = f.readline() # Discard header
                # Loop through data blocks
                for i in range(15): # Density points
                    buff = []    # Buffer
                    for j in range(12): # Temperature points
                        l = f.readline().split() 
                        # Arrange temperature points into 1D array
                        if j not in [0, 11]:
                            for x in l:
                                if 'E' not in x: 
                                    x = x.replace('-', 'E-').replace('+', 'E+')
                                buff.append(float(x))
                    # Append temperature points to reactions['IONIZ/REC'] 
                    # for each density point
                    reactions[k].append(buff)
                reactions[k] = transpose(array(reactions[k])) 
                                        # Convert 2D list into array(T,n)


    def get_coeff(self, database, reaction):
        ''' Returns the cofficients of reaction in database
            
        Parameters
        ----------
        database : string
            database in which to look for reaction
        reaction : string
            reaction handle/ID to retrieve from database

        Returns
        -------
        Reaction object
        '''
        return(self.reactions[database][reaction])



                



