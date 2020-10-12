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
            self, amjuel='amjuel.tex', hydhel='hydhel.tex', 
            h2vibr='h2vibr.tex', ADAS='ich0-1.dat', UE='ehr1.dat', 
            path='.'):
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
        self.reactions= {
                'AMJUEL': { 'settings' : [500, ['b0','0','a0','h0','p0','k0'],
                            ['a0','h0','p0','k0'], 45], 'path' : amjuel },
                'HYDHEL': {'settings' : [150, ['b0','0','a0','h0'], 
                            ['a0','h0'], 80], 'path' : hydhel },
                'H2VIBR': {'settings' : [0, ['b0','0','a0','h0'], ['a0','h0'],
                            20], 'path' : h2vibr },
                'UE': {},
                'ADAS': {} }
        # For each data point, add the reactions to the appropriate dictionary
        for rate in ['AMJUEL', 'H2VIBR', 'HYDHEL']:
            if self.reactions[rate]['path'] is not None:
                try:
                    self.read_EIRENE(self.reactions[rate]['path'],
                            self.reactions[rate], 
                            self.reactions[rate]['settings'], path=path)
                except:
                    print(  'Database {} not found in {}/{}. Omitting'
                            ''.format(rate, path, 
                            self.reactions[rate]['path']))
        if ADAS is not None:
            try:
                self.read_ADAS(ADAS, self.reactions['ADAS'], path=path)
            except:
                print(  'Database ADAS not found in {}/{}. Omitting'
                        ''.format(path, ADAS))
        if UE is not None: 
            try:
                self.read_UE(UE,self.reactions['UE'],path=path)
            except:
                print(  'Database UE not found in {}/{}. Omitting'
                        ''.format(path, UE))


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
        from numpy import zeros
        from csv import reader
        # Script parsing Hydhel.tex
        lst = []
        book = []
        # Open reader object
        rd = reader(open('{}/{}'.format(path,fname), 'rt'), delimiter=' ')
        # Store book in list
        for row in rd:
            lst.append(row)
        # Strip all empty arrays
        lst = list(filter(None, lst))

        for row in lst:
            book.append(list(filter(None, row)))

        i = settings[0] # Read from top
        while True:
            # Loop through rows looking for reactions
            if book[i][0] == "Reaction":
                reaction = book[i][1] # Store reaction name
                i += 1
                # Loop through reaction looking for coeffs
                while True:
                    if book[i][0] not in settings[1]:
                        i += 1
                     #   break
                    # Break if wrong fits
                    elif book[i][0] in settings[2]:
                        break
                    # We are in a T-fit
                    elif book[i][0] == 'b0':
                        coeff = []
                        # Parse the next three lines
                        for j in range(3):
                            for k in range(3):
                                coeff.append(float(
                                        book[i + j][1 + k*2].replace(
                                        ',', '').replace('D', 'E')))
                        i += 1
                        reactions[reaction.upper()] = coeff
                        break
                    # Wea re in a (T,E)-fit
                    elif book[i][0] == '0':
                        coeff = zeros((9, 9))
                        # Parse the whole data block in one
                        for j in range(3):
                            for k in range(9):
                                for m in range(3):
                                    coeff[k, j*3 + m]=float(
                                            book[i + k + j*9 + j*2][
                                            m+1].replace('D', 'E'))
                        # Store the coefficients
                        # TODO: figure out better way to kill off ne,T fits??
                        if reaction not in ['2.2.14', '2.0l2']:
                            reactions[reaction.upper()] = coeff
                        
                        i += 9 + 3*3 + 2*2 # Set counter after block
                        break
                    if i >= len(book) - settings[3]: # Omit last lines 
                        break

            i += 1
            if i >= len(book) - settings[3]: # Omit last lines
                break


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



                



