# CRUMPET input parser class readinput.py
# Separated from CRUMPET.py by holm10
# Changelog
# 200205 - Separated from CRUMPET.py #holm10
from .crm import Crm
from .ratedata import RateData
from cycler import cycler
import matplotlib.pyplot as plt
# Colorblind-compliant color cycle
colors = [  '#377eb8', 
            '#ff7f00', 
            '#4daf4a',
            '#f781bf',
            '#a65628',
            '#984ea3',
            '#999999',
            '#e41a1c',
            '#dede00',
            ]*4

lines = (['-' for x in range(9)]
        +['--' for x in range(9)]
        +[':' for x in range(9)]
        +['-.' for x in range(9)])

plt.rc('axes', prop_cycle=(cycler(color=colors) +
                  cycler(linestyle=lines)))



class Crumpet(Crm, RateData):
    '''CRUMPET container object for collsional-radiative models

    This class is used to read and parse the input file, the parts
    which are passed onto the CRM class that sets up the CRM itself.
    Inherits the classes Crm and RateData.
    
    ...

    Attributes
    ----------
    ver : srt
        a string defining the code version
    '''

    # TODO: Store input file name and include it in all output files!

    def __init__(
            self, fname='input/CRUMPET.dat', path='./', vmax=14, nmax=8,
            verbose=False, NP=2, isotope='H', mass=1):
        '''
        Parameters
        ----------
        fname : str, optional (default: 'input/CRUMPET.dat')
            path to CRUMPET input file relative to parameter path
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
        isotope : string, optional (default: 'H')
            The handle for the isotope being used. Used to identify 
            atomic and molecular species: molecules are taken to be
            {isotope}2{groundstatemol} and atoms to be 
            {isotope}{groundstateatom}
        mass : float, optional (default: 1)
            The isotope mass in AMU, used for calculations in Crm 
        '''
        from numpy import zeros, pi
        from os import getcwd

        self.ver = 'V1.0'
        self.ev = 1.602e-19
        self.fname  = fname

        def build_reactions(reaction_id):
            if reaction_id:
                return {reaction_id[0]: build_reactions(reaction_id[1:])}
            return []


        def merge_reactions(dct, merge_dct):
            from collections import Mapping
            """ Adapted from https://gist.github.com/angstwad/bf22d1822c38a92ec0a9
            """
            for k, v in merge_dct.items():
                if (k in dct and isinstance(dct[k], dict)
                        and isinstance(merge_dct[k], Mapping)):
                    merge_reactions(dct[k], merge_dct[k])
                else:
                    dct[k] = merge_dct[k]

        def parse_file(path, fname, data):
            ''' Parses file to dict '''
            with open('{}{}'.format(path, fname)) as f:  # Open the file
                for l in f:
                    # Omit empty and comment lines
                    l=l.strip()
                    if len(l)>0 and l[0] != '#':
                        l = l.split('#')[0].strip() # Discard comments
                        if '**' in l[:3]: # Card
                            card = l.split()[1].upper()
                            if card not in list(data):
                                data[card] = {}
                            subcard = None
                        elif '*' in l[:3]: # Subcard
                            if l.split()[1].strip().upper() == 'INCLUDE':
                                fnamepath = ' '.join(l.split()[2:])
                                 
                                parse_file('', ' '.join(l.split()[2:]), data)
                            elif card.upper() == 'REACTIONS':
                                if isinstance(data[card],dict):
                                    data[card] = []
                                data[card].append([l[2:]])
                                subcard = True
                                '''
                                try:
                                    [subcard, h123, reaction] = l.split()[1:]
                                    EIRreaction = reaction
                                    merge_reactions(data[card], build_reactions(l.split()[1:]))
                                except:
                                    # Account for name not matching reaction 
                                    # EIRENE database
                                    [subcard, h123, reaction, 
                                            EIRreaction] = l.split()[1:]
                                    merge_reactions(data[card], build_reactions(l.split()[1:-1]))
                                print(card, subcard)
                                data[card][subcard][h123][reaction].append(EIRreaction)'''
        
                            else:
                                subcard = l.split()[1]
                                data[card][subcard] = [] 
                                if len(l.split()) > 2:
                                    data[card][subcard] = l.split()[2]
                        else:
                            if subcard is None: # Direct data input w/o subcard
                                data[card][l.split()[0].strip()] = \
                                        ' '.join(l.split()[1:])
                            else:
                                if card.upper() == 'REACTIONS':
                                    #data[card][subcard][h123][reaction].append(l)
                                    data[card][-1].append(l)
                                else:
                                    data[card][subcard].append(l)


        # Parse the input file into the dict data
        data = {} # Hierarchy: Card - Subcard - Data
        parse_file(path, fname, data)

        ''' Store required input parameters '''
        # Species
        species = {}
        for key, value in data.pop('SPECIES',None).items():
            species[key] = {'V':0} # Set default value
            if len(value) > 0:
                for val in value:
                    # Update V if set in input
                    if val.split()[0].upper() == 'V':
                        species[key]['V'] = float(val.split()[1].strip())
        
        # Reactions
        rlist = data.pop('REACTIONS', None)



        ''' Store optional input parameters '''
        # Settings
        settings = data.pop('SETTINGS',{})
        # Rates
        rates = {}
        if 'RATES' in list(data):
            rates = data.pop('RATES', None)
        RateData.__init__(self, rates, path)
                
        # Background
        if 'BACKGROUND' in list(data):
            bg = {}
            for key, value in data.pop('BACKGROUND',None).items():
                bg[key] = {'V':0} # Default potential level
                if len(value) > 0:
                    for val in value:
                        # Update V if set in input
                        if val.split()[0].upper() == 'V':
                            bg[key]['V'] = float(val.split()[1].strip())
        else:
            bg = {'e':{'V':0}, 'p':{'V':0}} # Default background

        for key, values in data:
            print('Unrecognized card {}.\nIgnoring.'.format(setting))


        # %%%% Set up the CRM %%%%
        Crm.__init__(self, species, bg, rlist, path, self.reactions, settings)
     




    def write_ue_rates( self, fname='crumpet', E=0.1, Tm=0, 
                        groundstate=['(n=1)','(v=0)']):
        ''' Writes UEDGE-compatible tab-delimited rate files

        This script creates rates for H2 dissociation and H creation due
        to H2. The units are s**-1: hence, the particle balance is 
        determined by multiplying the H2 rates by 2, and the difference
        is due to MAR/MAI reactions. The files are written to a 
        DEGAS-style tab-delimited files with 60 temperature points and 
        15 density points distributed in log-log space according to 
        T[i]=10**(-1.2+i/10) eV for i in [0,59] and 
        n[j]=10**(10+j/2) cm**-3 for j in [0,14]. The files are saved in 
        Crumpet.path/output as {fname}_(E/n)rates.dat. To be compatible
        with UEDGE only rates for Np=2, where atoms and molecules are 
        included are presently possible. 

        Parameters
        ----------
        fname : str, optional (default: 'crumpet')
            prefix to tabulated rate data files
        E : floar, optional (default: 0.1)
            assumed temperature of the molecules for rates it Tm=False
        Tm : float, optional (default: 0)
            local molecular temperature for calculating energy losses 
            associated with molecules. Defaults to E if Tm=0.
        groundstate : lst of str, optional (default=['(n=1)', '(v=0)'])
            handles for the atomic and molecular species ground state
            identifier, in the order [atom, molecule]. Used to identify 
            the reactions associated with the atomic and molecular 
            ground states, in combination with Crm.isotope.

        Returns
        -------
        None
        '''

        from numpy import zeros,sum,array
        from os import mkdir
        from datetime import datetime
        from getpass import getuser
        from tqdm import tqdm

        # Parse into full handles
        h0h2=['{}{}'.format(self.isotope, groundstate[0]),
                '{}2{}'.format(self.isotope, groundstate[1])]  
        # Ensure that the CRM format coincides with the assumed format
        if list(self.species)[:2] != h0h2:
            print('This model assumes the first species to be H0 and the '
                  'second to be H2 in their electronic and vibrational '
                  'ground states. The first input species do not match '
                  '{} and {}, respectively.\nAborting.'
                  ''.format(h0h2[0],h0h2[1]))
            return

        elif self.Np != 2:
            print('Presently, UEDGE only considers two neutral hydrogenic '
                  'species, H0 and H2. The P-space of the model is larger '
                  'than this and, thus, is incompatible with the UEDGE '
                  'model. No rates can be written for Np!=2.\nAborting.')
            return

        # Create rate coefficient matrices
        ret = zeros((15, 60, self.Np, self.Np))
        ext = zeros((15, 60, self.Np))
        retE = zeros((15, 60, 5, self.Np))
        extE = zeros((15, 60, 5))

        # Calculate the Greenland (Np) space rates for the T-n space
        print('Creating UEDGE rate output')
        for i in tqdm(range(15)): # 15 density points in log-log space
    
            for j in range(60): # 60 temperature points in log-log space

                Te = 10**(-1.2 + j/10)
                ne = 10**(10 + 0.5*i) # Convert to ev and cm**-3

                # Store density rate to matrix
                ret[i,j,:,:], ext[i,j,:], _ = self.gl_crm(
                        Te, ne)
                
                # Calculate energy rates
                U = self.Sgl(Te, ne, Te, ne, E, Tm, write=False)
                # Include T-losses or not?

                # Sum over species to get net change
                retE[i, j,:,:] = array(
                        [sum(U[0][0], axis=0), sum(U[1][0], axis=0),
                        sum(U[2][0], axis=0), sum(U[3][0], axis=0),
                        sum(U[4][0], axis=0)])

                # Same for external energy source/sink
                extE[i, j,:] = array(
                        [sum(U[0][1], axis=0), sum(U[1][1], axis=0),
                        sum(U[2][1], axis=0), sum(U[3][1], axis=0),
                        sum(U[4][1], axis=0) ])
        
        # Ensure output directory exists 
        try:
            mkdir('{}/output'.format(self.path))
        except:
            pass
        
        # Write UEDGE data
        print(' Writing UEDGE reaction rates to {}'.format(fname + '_nrates'))
        with open(
                '{}/output/{}.dat'.format(self.path,fname + '_nrates'),'w') as f:
            # Headers and positions in Greenland matrices
            setups= [   [   ' H0 depl. Rate(jt,jn) (s**-1)  '
                            'Te(jt) = 10**(-1.2 + (jt-1)/10) jt= 1,60\n', 
                            ret[:,:, 0, 0]    ],
                        [   ' H2->H0 Rate(jt,jn) (s**-1)\n', 
                            ret[:,:, 0, 1]    ],
                        [   ' H2 depl. Rate(jt,jn) (s**-1)\n', 
                            ret[:,:, 1, 1]    ],
                        [   ' H2 creation Rate(jt,jn) (s**-1)\n', 
                            ret[:,:,1,0]    ],
                        [   ' H0 external source Rate(jt,jn) (s**-1)\n',
                            ext[:,:,0]      ],
                        [   ' H2 external source Rate(jt,jn) (s**-1)\n', 
                            ext[:,:,1]      ]
                    ]
      
            # Loop through all species contributions
            for m in range(len(setups)): 
                print('    Writing rate {}/{}'.format(m+1, len(setups)))
                f.write(setups[m][0])
                for j in range(15): # Write each density block
                    out=''
                    f.write(' jn =   {}'.format(j + 1)\
                            +(j == 0)*'; jt = 1 -> 60 going by '
                            'rows   ne(jn) = 10**(10 + 0.5*(jn-1)) jn=1,15' + '\n')

                    for k in range(60): # Store temperature data to string  
                        out+='{:1.5E}'.format(setups[m][1][j,k]).rjust(13,' ')
                        # Split string into ten rows
                        if (k + 1)/6 in range(1, 11):  
                            out += '\n'
                    f.write(out + '\n')   # Write the data to file
            f.write('\n{} by user {} using CRUMPET {}'.format(
                    datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    getuser(), self.ver))
            f.write('\nCalculated from file "{}"'.format(self.fname))


        print(' Writing UEDGE energy rates to {}'.format(fname + '_Erates'))
        with open('{}/output/{}.dat'.format(
                self.path, fname + '_Erates'),'w') as f:

            # Headers and positions in Greenland matrices
            setups = [
                    [' e-loss (H) (eV s**-1)  Te(jt) = 10**(-1.2 + '
                    '(jt-1)/10) jt= 1,60\n', retE[:,:,0,0]],
                    [' e-loss (H2) (eV s**-1)\n', retE[:,:,0,1]],
                    [' ext. e-loss (eV s**-1)\n', extE[:,:,0]],
                    [' ia source (H) (eV s**-1)\n', retE[:,:,1,0]],
                    [' ia source (H2) (eV s**-1)\n', retE[:,:,1,1]],
                    [' ext ia source (eV s**-1)\n', extE[:,:,1]],
                    [' Epot source (H) (eV s**-1)\n', retE[:,:,2,0]],
                    [' Epot source (H2) (eV s**-1)\n', retE[:,:,2,1]],
                    [' ext Epot source (eV s**-1)\n', extE[:,:,2]],
                    [' rad source,a (H) (eV s**-1)\n', retE[:,:,3,0]],
                    [' rad source,a (H2) (eV s**-1)\n', retE[:,:,3,1]],
                    [' ext rad source,a (eV s**-1)\n', extE[:,:,3]],
                    [' rad source,m (H) (eV s**-1)\n', retE[:,:,4,0]],
                    [' rad source,m (H2) (eV s**-1)\n', retE[:,:,4,1]],
                    [' ext rad source,m (eV s**-1)\n', extE[:,:,4]],
                ]
            # Loop through all species contributions
            for m in range(len(setups)): 
                print('    Writing rate {}/{}'.format(m + 1,len(setups)))
                f.write(setups[m][0])
                for j in range(15): # Write each density block
                    out = ''
                    f.write(' jn =   {}'.format(j+1)+(j==0)*
                            '; jt = 1 -> 60 going by rows   ne(jn) '
                            '= 10**(10 + 0.5*(jn-1)) jn=1,15'+'\n')
                    for k in range(60): # Store temperature data to string  
                        out += '{:1.5E}'.format(setups[m][1][j,k]).rjust(
                                13, ' ')
                        if (k + 1)/6 in range(1, 11): # Split string into ten rows
                            out += '\n'
                    f.write(out + '\n')   # Write the data to file
            f.write('\n {} by user {} using CRUMPET {}'.format(
                    datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    getuser(), self.ver))
            f.write('\nCalculated from file "{}"'.format(self.fname))


    def plot_depletion(
            self, key, n, Te, Ti=None, color='r', ax=None, E=0.1,
            figsize=(10,6.25), plot='loglog', **kwargs):
        ''' Plots the depletion rate for a given species

        Plots a curve or a set of curves for the depletion rate of the 
        species specified by key. If both n and Te are lists, a set of
        T-dependent rate curves are plotted. If either n or Te is a 
        float, a plot of as a function of the list is created. The 
        abcissa is Te regardless whether Ti is defined or not

        Parameters
        ----------
        key : str
            handle for the species to plot depletion rate of
        n : ndarray/list of floats or float
            plasma density values in cm**-3 
        Te : ndarray/list of floats or float
            electron temperature values in eV
        Ti : ndarray/list of floats or float, optional (default: None)
            ion temperature values in eV. Must be equal in length to
            Te. If None, Ti=Te is assumed
        ax : axes handle, optional (default: None)
            axes to plot curves on. If None, a new figure is created
        figsize : tuple (float,float), optional (default: (10,6.25))
            sets figure's figsize
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        plot : str, optional (default: plot)
            plot layout specifier, see documentation of plotax.
        **kwargs
            passed to ax.plot
    
        Returns
        -------
        matplotlib.pyplot.figure object
            if ax=None, otherwise returns None
        '''
        from matplotlib.pyplot import figure
        from numpy import ndarray
        if ax is None:
            fig = figure(figsize=figsize)
            ax = fig.add_subplot(111)
        if (not isinstance(n, list)) and (not type(n) is ndarray):
            n = [n]
        if (not isinstance(Te, list)) and (not type(Te) is ndarray):
            Te = [Te]
        if Ti is None: 
            Ti = Te
        elif (not isinstance(Ti, list)) and (not type(Ti) is ndarray):
            Ti = [Ti]
        if len(Te) != len(Ti):
            print('Te and Ti of different length: Te and Ti must be of\n'
                  'equal length. Aborting.')
            return 
        if len(Te) != len(Ti):
            print('Te and Ti must be of equal length!\nAborting.')
            return
        if len(Te) == 1:
            for i in range(len(Te)):
                y = []
                for k in range(len(n)):
                    y.append(self.species_rate(
                            key, Te[i], n[k], Ti[i], n[k],E))
                self.pickplot(ax, plot)(n, y, **kwargs)
        else:
            for k in range(len(n)):
                y = []
                for i in range(len(Te)):
                    y.append(self.species_rate(
                            key, Te[i], n[k], Ti[i], n[k], E))
                self.pickplot(ax, plot)(Te, y, **kwargs)
        try:
            fig.show()
            return fig 
        except:
            pass


    def plot_crm_dt(
            self, t, Te, ne, E=0.1, Ti=None, ni=None, Srec=True, Qres=True, 
            Tm=0, gl=False, density=True, fig=None, n0=None, title=None,
            figsize=(7+7,7*1.618033-3), savename=None, figtype='png',plotspecies=None,
            conservation=False, plot='plot', ext=None, plottot=False, **kwargs):
        ''' Plots the CRM results for neutrals in a static background plasma

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
        density : boolean, optional (default: False)
            Only returns density rates if True, otherwise evaluates
            both density and energy rates
        fig : matplotlib.pyplot.figure, optional (default: None)
            Figure which to plot on. By default creates a new figure and 
            returns the figure handle
        figsize : tuple (float, float), optional (default: (14, 8.33))
            Size of the created figure
        savename : string, optional (default: None)
            Name to save figure by. Figures are saved to 
            {Crumpet.path}/output{savename}.{figtype}. By default, 
            the figure is not saved.
        figtype : string, optional (default: 'png')
            filetype to save figure as. See matplotlib.pyplot.savefig
            for options
        title : string, optional (default: None)
            Figure title
        conservation : boolean, optional (default: False)
            Switch to also plot particle and energy conservation. 
            Particle conservation does not account for 
            recombination/ionization gains/losses
        ext : list of floats, optional (default: None)
            External sinks and sources

        Returns
        -------
        matplotlib.pyplot.figure
            if fig=None, else returns None
        '''
        # TODO: revise coding
        from matplotlib.pyplot import figure
        from numpy import sum, zeros, linspace, split, mod, concatenate, logspace, log10, array
        from numpy.linalg import inv
        from os import mkdir
        # TODO: Allow manual definition of ylim       
        # Get axis handle or create figure
        try:
            ax = fig.get_axes()[0]
        except:
            fig = figure(figsize=figsize)
        # If no initial distribution is requested, use the one specified in the input
        if n0 is None: 
            n0 = self.n0()
        if (Qres is True) and (gl is True): 
            print(  'Cannot resolve Q-species when Greenland-approximation is '
                    'used.\nAborting.')
            return
        # Solve the problem up to T
        # TODO: fix addext
        ft = self.solve_crm(t, Te, ne, Ti, ni, E, Tm, Srec, gl=gl, 
                            n=n0, Qres=Qres, densonly=density, **kwargs)#, addext=ext) 
        # TODO: allow passing of additional sources
        if ('logx' in plot) or ('loglog' in plot):
            t = logspace(-8, log10(t), 250)
        else:
            t = linspace(0, t, 250)
        xt = t*1e3
        # Solve densities only
        if density is True:
            if plotspecies is None:
                nt = ft.sol(t)
                legend = self.slist
            else:
                inds = []
                for s in plotspecies:
                    inds.append(self.slist.index(s))
                nt = ft.sol(t)[inds]
                legend = list(array(self.slist)[inds])
            fig.subplots_adjust(hspace=0, right=.8 + 0.13*Qres)
            ax = fig.add_subplot(111 + 100*conservation + 100*Qres)
            ylabel = r'Particles [$10^{{{}}}$]'
            self.plotax(ax, xt, nt, ylabel=ylabel, plottot=plottot,
                    plotmax=self.Np*(not Qres) + 1e3*Qres, plot=plot, legend=legend)
            if plottot:
                self.plottotpart(ax, n0)
            if Qres:
                ax.legend(  ncol=9 - 2*(not conservation), frameon=False,
                            bbox_to_anchor=(1.0, -0.4),
                            prop={'size':10 + 2*(not conservation)})
            else:
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                if conservation is True:
                    ax.set_xticklabels([])
            if conservation:
                ax = fig.add_subplot(212 + 101*Qres)
                ylabel = r'Particle balance [$10^{{{}}}$]'
                self.plotparterror(ax, xt, nt, ylabel=ylabel, color='k')
        # Solve densities and energy
        else:
            N = len(self.species)
            [Ee, Eia, Epot, Erada, Eradm, nt] = split(ft.sol(t),
                    [int(x) for x in linspace(N**Qres, 5*N**Qres, 5)], axis=0)
            if Qres is True:
                fig.subplots_adjust(hspace=0)
                data=[nt, -Ee, Eia, Epot, Erada, Eradm]
                ylabel=[
                        r'Particles [$10^{{{}}}$]',
                        r'$\rm{{ S_{{e,sink}} }}$ [$10^{{{}}}$ J]',
                        r'$\rm{{ S_{{i/a}} }}$ [$10^{{{}}}$ J]',
                        r'$\rm{{ S_{{pot}} }}$ [$10^{{{}}}$ J]',
                        r'$\rm{{ S_{{rad,a}} }}$ [$10^{{{}}}$ J]',
                        r'$\rm{{ S_{{rad,m}} }}$ [$10^{{{}}}$ J]']
                for i in range(len(data)):
                    ax = fig.add_subplot(400 + 100*conservation + 2*10 + i + 1)
                    self.plotax(ax, xt, data[i], nuclei=False, 
                            ylabel=ylabel[i], plot=plot)
                    if i < 4:
                        ax.set_xticklabels([])


                ax.legend(  ncol=11 - 4*(not conservation), frameon=False,
                            bbox_to_anchor=(1.03, -0.4),
                            prop={'size':8 + 4*(not conservation)})
                # Plot large grid with individual species
                if conservation:
                    ax = fig.add_subplot(515)
                    ylabel = r'Energy balance [$10^{{{}}}$ J]'
                    self.plotEerror(ax, xt,
                            self.ev*concatenate((Ee, Eia, Epot, Erada, Eradm)),
                            ylabel=ylabel)
                    ax.plot([], [], 'k-')
                    self.plotparterror(ax, xt, nt, right=True, color='k',
                            ylabel=r'Neutral nuclei conservation '
                            '[$10^{{{}}}$]',)
                    ax.legend(['Energy', 'Neutral nuclei'],
                            title='Conservation')
                self.plottotpart(fig.get_axes()[0], n0)
                self.setgrid(fig)
            else:
                fig.subplots_adjust(hspace=0, right=0.8)
                ax1 = fig.add_subplot(211 + 100*conservation)
                ylabel = r'Particles [$10^{{{}}}$]'
                self.plotax(ax1, xt, nt, ylabel=ylabel, plotmax=self.Np, 
                            plot=plot)
                self.plottotpart(ax1, n0)
                ax1.set_xticklabels([])
                ax1.legend( title='Particle balance', bbox_to_anchor=(1.05, 1),
                            loc='upper left')

                ax2 = fig.add_subplot(212 + 100*conservation)
                ylabel = r'Energy [$10^{{{}}}$ J]'
                self.plotax(ax2, xt,
                        self.ev*concatenate((-Ee, Eia, Epot, Erada, Eradm)),
                        plottot=False, ylabel=ylabel, plot=plot)
                ax2.legend([r'$\rm S_{e,loss}$',r'$\rm S_{i/a}$',
                            r'$\rm S_{pot}$',r'$\rm S_{rad,a}$',
                            r'$\rm S_{rad,m}$'],
                            title='Energy sinks/sources',
                            bbox_to_anchor=(1.05, 1), loc='upper left')
                if conservation is True: 
                    ax2.set_xticklabels([])
                if conservation:
                    ax3 = fig.add_subplot(213 + 100*conservation)
                    self.plotEerror(ax3, xt,
                            self.ev*concatenate((Ee, Eia, Epot, Erada, Eradm)),
                            ylabel=r'Energy conservation [$10^{{{}}}$ J]')
                    ax3.plot([], [], 'k-')
                    self.plotparterror(ax3, xt, nt, right=True, color='k',
                            ylabel=r'Neutral nuclei conservation '
                            '[$10^{{{}}}$]')
                    ax3.legend([    'Energy','Neutral nuclei'],
                                    title='Conservation',loc='upper left',
                                    bbox_to_anchor=(1.05, 1)    )
                self.setgrid(fig)
        if title is not None:
            fig.suptitle(title)
        # Save to casepath/output/figs/figname.type if figname set
        if savename is not None:
            # Make sure that the direcory exists
            try:
                mkdir('{}/output/figs'.format(self.path))
            except:
                pass
            # Save
            fig.savefig('{}/output/figs/{}.{}'.format(
                        self.path,savename,figtype), dpi=300, edgecolor=None,
                        format=figtype, bbox_inches='tight')
        try:
            fig.show() # Show fig
            return fig
        except: 
            return


    def plot_species_dt(
            self, t, Te, ne, species, E=0.1, Ti=None, ni=None, Srec=True, Qres=True, 
            Tm=0, gl=False, density=False, fig=None, n0=None, title=None, ylim=(None,None),
            figsize=(7+7,7*1.618033-3), savename=None, figtype='png', plot='semilogx', legend=None,
            xlabel=None):
        ''' 
        TODO: write documentation
        Plots the CRM results for neutrals in a static background plasma

        Parameters
        ----------
        t : float
            end-time of simulation in s
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        species : list of strings
            list containing species handles to be plotted
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
        density : boolean, optional (default: False)
            Only returns density rates if True, otherwise evaluates
            both density and energy rates
        fig : matplotlib.pyplot.figure, optional (default: None)
            Figure which to plot on. By default creates a new figure and 
            returns the figure handle
        figsize : tuple (float, float), optional (default: (14, 8.33))
            Size of the created figure
        savename : string, optional (default: None)
            Name to save figure by. Figures are saved to 
            {Crumpet.path}/output{savename}.{figtype}. By default, 
            the figure is not saved.
        figtype : string, optional (default: 'png')
            filetype to save figure as. See matplotlib.pyplot.savefig
            for options
        title : string, optional (default: None)
            Figure title

        Returns
        -------
        matplotlib.pyplot.figure
            if fig=None, else returns None
        '''
        from matplotlib.pyplot import figure
        from numpy import sum, zeros, linspace, split, mod, concatenate, array
        from numpy.linalg import inv
        from os import mkdir
        # TODO: Allow manual definition of ylim       
        # Get axis handle or create figure

        try:
            ax = fig.get_axes()[0]
            fig.subplots_adjust(bottom=.3)
            prefig = True
        except:
            prefig = False
            fig = figure(figsize=figsize)
            fig.subplots_adjust(hspace=0, right=.8 + 0.13*Qres)
            ax = fig.add_subplot(111 + 100*Qres)
        # If no initial distribution is requested, use the one specified in the input
        if n0 is None: 
            n0 = self.n0()
        if (Qres is True) and (gl is True): 
            print(  'Cannot resolve Q-species when Greenland-approximation is '
                    'used.\nAborting.')
            return
        # Solve the problem up to T
        ft = self.solve_crm(t, Te, ne, Ti, ni, E, Tm, Srec, gl=gl, 
                            n=n0, Qres=Qres, densonly=density) 
        inds = []
        for sp in species:
            inds.append(self.get_species_ind(sp))

        # Solve densities only
        if density is True:
            nt = [] 
            for i in inds:
                nt.append(ft.y[i])
            nt=array(nt)
            ylabel = r'Particles [$10^{{{}}}$]'
            self.plotax(ax, ft.t*1e3, nt, ylabel=ylabel, xlabel=xlabel, plottot=False,
                    plotmax=self.Np*(not Qres) + 1e3*Qres, plot=plot,
                    xlim=(1e-8,1e2), ylim=ylim, legend=legend)
            if Qres:
                ax.legend(  ncol=7, frameon=False,
                            bbox_to_anchor=(1.0, -0.4+0.25*prefig),
                            prop={'size':12})
            else:
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        '''
        # Solve densities and energy
        else:
            N = len(self.species)
            [Ee, Eia, Epot, Erada, Eradm, nt] = split(ft.sol(t),
                    [int(x) for x in linspace(N**Qres, 5*N**Qres, 5)], axis=0)
            if Qres is True:
                fig.subplots_adjust(hspace=0)
                data=[nt, -Ee, Eia, Epot, Erada, Eradm]
                ylabel=[
                        r'Particles [$10^{{{}}}$]',
                        r'$\rm{{ S_{{e,sink}} }}$ [$10^{{{}}}$ J]',
                        r'$\rm{{ S_{{i/a}} }}$ [$10^{{{}}}$ J]',
                        r'$\rm{{ S_{{pot}} }}$ [$10^{{{}}}$ J]',
                        r'$\rm{{ S_{{rad,a}} }}$ [$10^{{{}}}$ J]',
                        r'$\rm{{ S_{{rad,m}} }}$ [$10^{{{}}}$ J]']
                for i in range(len(data)):
                    ax = fig.add_subplot(400 + 2*10 + i + 1)
                    self.plotax(ax, xt, data[i], nuclei=False, 
                            ylabel=ylabel[i], plot=plot, plottot=False,
                            xlim=(1e-8,t),ylim=ylim)
                    if i < 4:
                        ax.set_xticklabels([])


                ax.legend(  ncol=7, frameon=False,
                            bbox_to_anchor=(1.03, -0.4),
                            prop={'size':12})
                # Plot large grid with individual species
                self.setgrid(fig)
            else:
                fig.subplots_adjust(hspace=0, right=0.8)
                ax1 = fig.add_subplot(211)
                ylabel = r'Particles [$10^{{{}}}$]'
                self.plotax(ax1, xt, nt, ylabel=ylabel, plotmax=self.Np, plot=plot,
                        plottot=False, xlim=(1e-8,1e2))
                ax1.set_xticklabels([])
                ax1.legend( title='Particle balance', bbox_to_anchor=(1.05, 1),
                            loc='upper left')

                ax2 = fig.add_subplot(212)
                ylabel = r'Energy [$10^{{{}}}$ J]'
                self.plotax(ax2, xt,
                        self.ev*concatenate((-Ee, Eia, Epot, Erada, Eradm)),
                        plottot=False, ylabel=ylabel, plot=plot, xlim=(1e-8,1e2))
                ax2.legend([r'$\rm S_{e,loss}$',r'$\rm S_{i/a}$',
                            r'$\rm S_{pot}$',r'$\rm S_{rad,a}$',
                            r'$\rm S_{rad,m}$'],
                            title='Energy sinks/sources',
                            bbox_to_anchor=(1.05, 1), loc='upper left')
                self.setgrid(fig)
        '''
        if title is not None:
            fig.suptitle(title)
        # Save to casepath/output/figs/figname.type if figname set
        if savename is not None:
            # Make sure that the direcory exists
            try:
                mkdir('{}/output/figs'.format(self.path))
            except:
                pass
            # Save
            fig.savefig('{}/output/figs/{}.{}'.format(
                        self.path,savename,figtype), dpi=300, edgecolor=None,
                        format=figtype, bbox_inches='tight')
        try:
            fig.show() # Show fig
            return fig
        except: 
            return


    def read_ue_crumpet(self, nratefile, eratefile):
        ''' **INCOMPLETE** Reads a CRUMPET UEDGE rate file and returns the rates
        
        Parameters
        ----------
        nratefile : string
            path to Crumpet density rate file relative to Crumpet.path
        Eratefile : string
            path to Crumpet energy rate file relative to Crumpet.path
        
        Returns
        -------
        dict
            the dict contains the atomic and molecular density rates
            and the energy rates (internal and external)

        '''
        from CRUMPET.reactions import Reaction
        rates = {}
        reactions = {}
        # Read the data from ratefile into dict
        datalist = ['H0_depl', 'H0_create', 'H2_depl', 'H2_create', 'H0_ext',
                'H2_ext']

        self.read_UE(nratefile, rates, datalist=datalist, path=self.path)
        # Create custom reaction
        reactions = {}
        print(rates.keys())
        for r in rates.keys():
#            reactions[r] = Reaction(r, rates[r], 'UE', {'reactants':['H2(v=0)'], 
#                    'fragments':['H2(v=0)'], 'K':'0'}, self.bg,
#                    self.species, self.isotope, self.mass)
            print(rates[r])
            reactions[r] =  Reaction('UE', 'UE', r, reactiondata, coeffs, self.bg, 
                self.species, self.isotope, self.mass, scale=scale)

        # Read the data from ratefile into dict
        datalist=['Hel', 'H2el', 'extel', 'Hia', 'H2ia', 'extia', 'HV', 'H2V',
                'extV', 'Hga', 'H2ga', 'extga', 'Hgm', 'H2gm', 'extgm',]
        try:
            self.read_UE(eratefile, rates, datalist=datalist, path=self.path)
            # Create custom reaction
            for r in rates.keys():
                rates[r] = Reaction(r, rates[r], 'UE', {'reactants':['H2(v=0)'], 
                                    'fragments':['H2(v=0)'], 'K':'0'}, self.bg,
                                    self.species, self.isotope, self.mass)
        except:
            pass
        return rates


    def plot_ue_nrate(
            self, ratefile, Te, ne, fig=None, figsize=(10,10/1.618033),
            linestyle='-', ylim=(1e-1,1e7), linewidth=3, labelapp='',
            savename=None, title=None, pretitle='', figtype='png', 
            color='k', plot='loglog', idx=0, xlim=None, ncol=3, fac=1,
            E=0.1, Tm=0):
        ''' **INCOMPLETE** Plots the density rate coefficients read from a crumpet rate file
        
        Parameters
        ----------
        ratefile : string
            path to Crumpet density rate file relative to Crumpet.path
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        Tm : float, optional (default: 0)
            local molecular temperature for calculating energy losses 
            associated with molecules. Defaults to E if Tm=0.
        fig : matplotlib.pyplot.figure, optional (default: None)
            figure which to plot on. By default creates a new figure and 
            returns the figure handle
        figsize : tuple (float, float), optional (default: (14, 8.33))
            size of the created figure
        linestyle : string, optional (default: '-')
            linestyle option for matplotlig.pyplot used in the plot
        savename : string, optional (default: None)
            Name to save figure by. Figures are saved to 
            {Crumpet.path}/output{savename}.{figtype}. By default, 
            the figure is not saved.
        ylim : len-2 tuple of floats, optional (default: (1e-7, 1e7))
            ordinate limits in cm**3 s**-1
        xlim : len-2 tuple of floats, optional (default: None)
            abscissae limits in cm**3 s**-1, defaults to default limits
        linewidth : float, optional (default: 3)
            width of lines plotted 
        labelapp : string, optional (default: '')
            appendix for plot labels
        title : string, optional 
        figtype : string, optional (default: 'png')
            filetype to save figure as. See matplotlib.pyplot.savefig
            for options
        title : string, optional (default: None)
            figure title. Defaults to temperature and density
        pretitle : string, optional (default: '')
            prefix to default title, not used if title is not None
        color : string, optional (default: 'k')
            color to be used in plot
        plot : sting, optional (default: 'loglog')
            string of the matplotlib.pyplot.Axes plot option to be used
        idx : int, optional (default: 0)
            Index of appearance of the data block in ratefile:
            0 - atom sink
            1 - atom source 
            2 - molecular sink
            3 - molecular source
            4 - external (recombination) atom source    
            5 - external (re-association) molecular source
        ncol : int, optional (default: 3)
            ncol keyword of legend: how many columns in legend
        fac : float, optional (default: 1)
            scaling factor for ordinate 

        Returns
        -------
        matplotlib.pyploy.Figure object
        '''
        from matplotlib.pyplot import figure
        from numpy import log10, array, ndarray
        from os import mkdir
        from CRUMPET.reactions import Reaction

        try:
            rates = {}
            # Read the data from ratefile into dict
            datalist = ['H0_depl', 'H0_create', 'H2_depl', 'H2_create', 
                    'H0_ext', 'H2_ext']
            self.read_UE(ratefile, rates, datalist=datalist)
            # Create custom reaction
            reactions = {}
            for r in rates.keys():
                reactions[r] = Reaction(r, '', rates[r], 'UE',
                        ['', 0, '', 0, [0, 0, 0, 0]])
        except:
            reactions=self.calc_ue_nrate(Te, ne, E, Tm)
        typ = ''
        lab = ''
        if idx in [0, 1, 4]:
            typ += r'Atomic hydrogen '
            lab += 'Atom. H '
        else:
            typ += r'Molecular hydrogen '
            lab += 'Mol. H '
        if idx in [0, 2]:
            typ += 'depletion'
            lab += 'sink'
        elif idx in [1, 3]:
            typ += 'creation'
            lab += 'source'
        else:
            typ += 'external source (EIR)'
            lab += 'ext. source'
        # Setup plot type
        if isinstance(Te, list) or isinstance(Te, ndarray):
            x = Te
            try:
                y = array([reactions[datalist[idx]].rate(
                        T, None, ne=ne) for T in Te])/(ne**(idx > 3))
            except:
                y = reactions[datalist[idx]]/(ne**idx > 3)
            xlabel = (r'Electron temperature [eV]')
            ex = int(log10(ne))
            mu = ne/(10**ex) 
            autotitle = pretitle+ r' {}, $\rm{{n_e}}$={}$\times10^{{{}}}$ '\
                                '$\rm{{cm^{{-3}}}}$'.format(typ, mu, ex)
        elif isinstance(ne, list): 
            x = ne
            y = array([reactions[datalist[idx]].rate(
                    Te, None, ne=n) for n in ne])
            if idx > 3:
                y = y/ne
            xlabel = (r'Electron density [$\rm{cm^{-3}}$]')
            autotitle = pretitle + r' {}, $\rm{{T_e}}$={} [eV]'.format(typ, Te)
        else:
            print('Reaction rate for process "{}" at Te={} and ne={}'.format(
                    datalist[idx], Te, ne))
            return reactions[datalist[idx]].rate(Te, None, ne)
        try:
            ax = fig.get_axes()[0]
        except:
            fig = figure(figsize=figsize)
            ax = fig.add_subplot(111)

        if plot == 'loglog':
            pl = ax.loglog
        elif plot == 'semilogy':
            pl = ax.semilogy
        elif plot == 'semilogx':
            pl = ax.semilogx
        elif plot == 'plot':
            pl = ax.plot
        else:
            print('Unrecognized plot "{}" requested! Aborting.'.format(plot))
        pl(x, fac*abs(y), linewidth=linewidth, color=color, 
                label=lab + labelapp, linestyle=linestyle)
        if title is None:
            ax.set_title(autotitle)
        else:
            ax.set_title(title)
        ax.legend(ncol=ncol, loc='best')
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
            fig.savefig('output/figs/{}.{}'.format(savename, figtype), dpi=300,
                    edgecolor=None, format=figtype, bbox_inches='tight')
        fig.show()
        return fig

    def calc_ue_nrate(self, Te, ne, E, Tm):
        ''' **INCOMPLETE** Calculates the Greenland CRM particle rate coefficients

        Includes external sources 

        Parameters
        ----------
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        Tm : float, optional (default: 0)
            local molecular temperature for calculating energy losses 
            associated with molecules. Defaults to E if Tm=0.

        Returns
        -------
        dict
            Returns a dict of rate coefficients, with the corresponding
            handles:
            H0_depl - atom sink
            H0_create - atom source
            H2_depl - molecular sink
            H2_create - molecular source
            H0_ext - external (recombination) atom source
            H2_ext - external (re-association) molecule sink
        '''
        from numpy import zeros, sum, array

        datalist = ['H0_depl', 'H0_create', 'H2_depl', 'H2_create', 'H0_ext',
                'H2_ext']
        ret = zeros((len(Te), self.Np, self.Np))
        ext = zeros((len(Te), self.Np))
        for i in range(len(Te)):
            T = Te[i]
            # Store to matrix
            ret[i,:,:], ext[i,:], _ = self.gl_crm(T, ne)
        setups = [
                ['H0_depl', ret[:,0,0]],
                ['H0_create', ret[:,0,1]],
                ['H2_depl', ret[:,1,1]],
                ['H2_create', ret[:,1,0]],
                ['H0_ext', ext[:,0]],
                ['H2_ext', ext[:,1]],
                ]
        reactions = {}
        for m in setups: # Loop through all species contributions
            reactions[m[0]] = m[1]
        return reactions


    def calc_ue_erate(self,Te,ne,E,Tm):
        ''' **INCOMPLETE** Calculates the Greenland CRM energy rate coefficients

        Includes external sources 

        Parameters
        ----------
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        Tm : float, optional (default: 0)
            local molecular temperature for calculating energy losses 
            associated with molecules. Defaults to E if Tm=0.

        Returns
        -------
        dict
            Returns a dict of rate coefficients, with the corresponding
            handles:
            Hel - electron energy sink/source due to atoms 
            H2el - electron energy sink/source due to molecules
            extel - external electron energy sink/source
            Hia - ion/atom energy sink/source due to atoms
            H2ia - ion/atom energy sink/source due to molecules
            extia - external ion/atom energy sink/source
            HV - potential energy sink/source due to atoms
            H2V - potential energy sink/source due to molecules
            extV - external potential energy sink/source    
            Hga - atomic line radiation power due to atoms
            H2ga - atomic line radiation power due to molecules
            extga - external atomic line radiation power 
            Hgm - molecular line radiation power due to atoms
            H2gm - molecular line radiation power due to molecules
            extgm - external molecular line radiation power 
        '''
        from numpy import zeros, sum, array

        datalist = ['Hel', 'H2el', 'extel', 'Hia', 'H2ia', 'extia', 'HV',
                    'H2V', 'extV', 'Hga', 'H2ga', 'extga', 'Hgm', 'H2gm',
                    'extgm',]
        retE = zeros((len(Te), 5, self.Np))
        extE = zeros((len(Te), 5))
        for i in range(len(Te)):
            T = Te[i]
            # Don't include erl1/erl2 radiation in the rates as these are handled by UEDGE
            U = self.Sgl(T, ne, T, ne, E, Tm, write=False) 
            # Include T-losses or not?
            retE[i,:,:] = array([
                    sum(U[0][0], axis=0), sum(U[1][0], axis=0), 
                    sum(U[2][0], axis=0), sum(U[3][0], axis=0), 
                    sum(U[4][0], axis=0)])
            extE[i,:] = array([
                    sum(U[0][1], axis=0), sum(U[1][1], axis=0), 
                    sum(U[2][1], axis=0), sum(U[3][1], axis=0), 
                    sum(U[4][1], axis=0) ])
        setups = [
                ['Hel', retE[:,0,0]],
                ['H2el', retE[:,0,1]],
                ['extel', extE[:,0]],
                ['Hia', retE[:,1,0]],
                ['H2ia', retE[:,1,1]],
                ['extia', extE[:,1]],
                ['HV', retE[:,2,0]],
                ['H2V', retE[:,2,1]],
                ['extV', extE[:,2]],
                ['Hga', retE[:,3,0]],
                ['H2ga', retE[:,3,1]],
                ['extga', extE[:,3]],
                ['Hgm', retE[:,4,0]],
                ['H2gm', retE[:,4,1]],
                ['extgm', extE[:,4]],
                ]
        reactions = {}
        for m in setups: # Loop through all species contributions
            reactions[m[0]] = m[1]
        return reactions


    def plot_ue_erate(
            self, ratefile, Te, ne, fig=None, origUE=False, idx=0, E=0.1, 
            Tm=0, figsize=(12,13/1.618033), linestyle='-', 
            ylim=(1e-17,1e-11), linewidth=2, labelapp='', savename=None,
            title=None, pretitle='', figtype='png', plot='loglog', ncol=3,
            xlim=None, colors=['darkcyan','b','m','r'], eion=5, ediss=10):
        ''' **INCOMPLETE** Plots the energy rate coefficients read from a crumpet rate file
        
        Plots the energy rate coefficients for the process with idx, as 
        explained below. Reads a CRUMPET energy rate file if origUE is 
        False, or emulates the original UEDGE calculation based on a 
        CRUMPET density rate file if origUE is True

        Parameters
        ----------
        ratefile : string
            path to Crumpet density rate file relative to Crumpet.path. 
            Energy file if origUE=False, density file if origUE=True
        Te : float
            electron temperature in eV
        ne : float
            electron density in cm**-3
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        Tm : float, optional (default: 0)
            local molecular temperature for calculating energy losses 
            associated with molecules. Defaults to E if Tm=0.
        fig : matplotlib.pyplot.figure, optional (default: None)
            figure which to plot on. By default creates a new figure and 
            returns the figure handle
        figsize : tuple (float, float), optional (default: (14, 8.33))
            size of the created figure
        linestyle : string, optional (default: '-')
            linestyle option for matplotlig.pyplot used in the plot
        savename : string, optional (default: None)
            Name to save figure by. Figures are saved to 
            {Crumpet.path}/output{savename}.{figtype}. By default, 
            the figure is not saved.
        ylim : len-2 tuple of floats, optional (default: (1e-7, 1e7))
            ordinate limits in cm**3 s**-1
        xlim : len-2 tuple of floats, optional (default: None)
            abscissae limits in cm**3 s**-1, defaults to default limits
        linewidth : float, optional (default: 3)
            width of lines plotted 
        labelapp : string, optional (default: '')
            appendix for plot labels
        title : string, optional 
        figtype : string, optional (default: 'png')
            filetype to save figure as. See matplotlib.pyplot.savefig
            for options
        title : string, optional (default: None)
            figure title. Defaults to temperature and density
        pretitle : string, optional (default: '')
            prefix to default title, not used if title is not None
        color : string, optional (default: 'k')
            color to be used in plot
        plot : sting, optional (default: 'loglog')
            string of the matplotlib.pyplot.Axes plot option to be used
        origUE : boolean, optional (default: False)
            switch to use the original UEDGE dissociation rate (True) 
            or not (False). calculated based on eion and ediss
        eion : float, optional (default: 5)
            Value of bbb.eion used in UEDGE: ion birth energy
        ediss : float, optional (default: 10)
            Value of bbb.ediss used in UEDGE: molecular dissociation 
            energy
        colors : list of strings, optiona 
                            (default:['darkcyan','b','m','r'])
            List of colors to be used in the plots
        idx : int, optional (default: 0)
            Index of appearance of the data block in ratefile:
            0 - electron energy sink/source due to atoms 
            1 - electron energy sink/source due to molecules
            2 - external electron energy sink/source
            3 - ion/atom energy sink/source due to atoms
            4 - ion/atom energy sink/source due to molecules
            5 - external ion/atom energy sink/source
            6 - potential energy sink/source due to atoms
            7 - potential energy sink/source due to molecules
            8 - external potential energy sink/source    
            9 - atomic line radiation power due to atoms
            10 - atomic line radiation power due to molecules
            11 - external atomic line radiation power 
            12 - molecular line radiation power due to atoms
            13 - molecular line radiation power due to molecules
            14 - external molecular line radiation power 
        ncol : int, optional (default: 3)
            ncol keyword of legend: how many columns in legend
        fac : float, optional (default: 1)
            scaling factor for ordinate 

        Returns
        -------
        matplotlib.pyploy.Figure object
        '''
        from matplotlib.pyplot import figure
        from numpy import log10, array, zeros, sum, ndarray
        from os import mkdir
        from CRUMPET.reactions import Reaction
        ev = 1.602e-19
        rates = {}
        try:
            if origUE is True:
                # Read the data from ratefile into dict
                datalist = ['H0_depl', 'H0_create', 'H2_depl', 'H2_create',
                            'H0_ext', 'H2_ext']
                self.read_UE(ratefile, rates, datalist=datalist)
            else:
                # Read the data from ratefile into dict
                datalist = ['Hel', 'H2el', 'extel', 'Hia', 'H2ia', 'extia',
                            'HV', 'H2V', 'extV', 'Hga', 'H2ga', 'extga', 'Hgm',
                            'H2gm', 'extgm',]
                self.read_UE(ratefile, rates, datalist=datalist)
            # Create custom reaction
            reactions = {}
            for r in rates.keys():
                reactions[r] = Reaction(r, '', rates[r], 'UE',
                        ['', 0, '', 0, [0, 0, 0, 0]])
        except:
            reactions = self.calc_ue_erate(Te, ne, E, Tm)
        title = ['Electron loss', 'Ion/atom source', 'Potential source',
                'Radiation source']
        lab = ['H0', 'H2', 'ext', 'H2,a', 'H2,m']
        try:
            ax = fig.get_axes()[0]
            oldax = True
        except:
            fig = figure(figsize=figsize)
            oldax = False
        for i in range(4):
            ax = fig.add_subplot(2, 2, i + 1)

            if plot == 'loglog':
                pl = ax.loglog
            elif plot == 'semilogy':
                pl = ax.semilogy
            elif plot == 'semilogx':
                pl = ax.semilogx
            elif plot == 'plot':
                pl = ax.plot
            else:
                print('Unrecognized plot "{}" requested! '
                        'Aborting.'.format(plot))
            if origUE:
                fac = (ev*(ediss*(i == 0) + 2*eion*(i == 1) 
                        + (ediss - 2*eion)*(i == 3)*(ediss > 2*eion)
                        - (ediss - 2*eion)*(i == 2)*(ediss < 2*eion)))
                # Setup plot type
                if isinstance(Te, list):
                    x = Te
                    y = array([reactions['H2_depl'].rate(
                            T, None, ne=ne)*fac for T in Te])
                    xlabel = (r'Electron temperature [eV]')
                    ex = int(log10(ne))
                    mu = ne/(10**ex) 
                    autotitle = title[i]
                    suptitle = pretitle + r'$\rm{{n_e}}$={}$\times10^{{{}}}$'\
                            ' $\rm{{cm^{{-3}}}}$'.format(mu, ex)
                elif isinstance(ne, list): 
                    x = ne
                    y = array([reactions['H2_depl'].rate(
                                        Te, None, ne=n)*fac for n in ne])
                    if idx > 3:
                        y = y/ne
                    xlabel = (r'Electron density [$\rm{cm^{-3}}$]')
                    autotitle = title[i]
                    supttitle = (pretitle + r' {}, $\rm{{T_e}}$={} ' 
                            + '[eV]'.format(typ, Te))
                else:
                    print('Reaction rate for process "{}" at Te={}'
                            ' and ne={}'.format(datalist[idx], Te, ne))
                    return reactions[datalist[idx]].rate(Te, None, ne)
                   
                pl(x, ((-1)**(i < 2))*y, linewidth=linewidth, color='k',
                        label='Orig. UE' + labelapp, linestyle=linestyle)
                ax.set_title(autotitle)
            else:
                if i != 3:
                    ii = 3
                else:
                    ii = 6
                for j in range(ii):
                    if ii == 6:
                        if j > 2:
                            jj = 3
                        else:
                            jj = j
                    else:
                        jj = j
                    # Setup plot type
                    if isinstance(Te, list) or isinstance(Te, ndarray):
                        x = Te
                        try:
                            y = array([reactions[datalist[i*3 + j]].rate(
                                    T, None, ne=ne)*ev for T in Te])/\
                                    (ne**(j == 2))
                        except:
                            y = reactions[datalist[i*3 + j]]*ev
                        xlabel = (r'Electron temperature [eV]')
                        ex = int(log10(ne))
                        mu = ne/(10**ex) 
                        autotitle = title[i]
                        suptitle = (pretitle+r'$\rm{{n_e}}$={}$\times10^{{{}}}$'
                                + ' $\rm{{cm^{{-3}}}}$'.format(mu, ex))
                    elif isinstance(ne, list): 
                        x = ne
                        y = array([reactions[datalist[i*3 + j]].rate(
                                Te, None, ne=n)*ev for n in ne])
                        if j == 2:
                            y = y/ne
                        xlabel = (r'Electron density [$\rm{cm^{-3}}$]')
                        autotitle = title[i]
                        supttitle = (pretitle + r' {}, $\rm{{T_e}}$={} '
                                + '[eV]'.format(typ, Te))
                    else:
                        print(  'Reaction rate for process "{}" at Te={}'
                                ' and ne={}'.format(datalist[idx], Te, ne))
                        return reactions[datalist[idx]].rate(Te, None, ne)
                    if j == 0:
                        tot = zeros((len(x),))
                    tot = tot + y
                    pl(x, ((-1)**(i == 0))*y, color=colors[jj], label=lab[jj]
                            +labelapp, linewidth=linewidth, 
                            linestyle=linestyle)
                    ax.set_title(autotitle)
                if (plot not in ['semilogy', 'loglog']) and (origUE is False):
                    pl(x, ((-1)**(i == 0))*tot, linewidth=linewidth, color='k',
                            label='Total' + labelapp, linestyle=linestyle)
            if oldax is False: 
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.07, box.width,
                         box.height * 0.95])
            if i == 3:
                ax.legend(loc='lower center', bbox_to_anchor=(-0.1, -0.35),
                        frameon=False, ncol=ncol)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(r'$\rm{S_E}$ [$\rm{W/(n_e\cdot V\quad)}$]')
            if xlim is not None:
                ax.set_xlim(xlim)
            ax.set_ylim(ylim)
#        fig.suptitle(suptitle)
        if savename is not None:
            try:
                mkdir('output/figs')
            except:
                pass
            fig.savefig('output/figs/{}.{}'.format(savename, figtype), dpi=300,
                    edgecolor=None, format=figtype, bbox_inches='tight')
        fig.show()
        return fig


    def spectrum(
            self, Te, ne, E=0.1, Ti=None, ni=None, Srec=True, n0=None, 
            intensity=None, units='l',write=False, norm=False, fig=None,
            figsize=(10,10/1.618033), xlim=None, linewidth=1, split=False,
            **kwargs):
        ''' Plots atomic and molecular spectra 

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
        fig : matplotlib.pyplot.figure, optional (default: None)
            figure which to plot on. By default creates a new figure and 
            returns the figure handle
        figsize : tuple (float, float), optional (default: (14, 8.33))
            size of the created figure
        units : string, optional (default: l)
            units of the plot ordinate
            ev - power in eV/s
            v - wavenumber cm**-1
            f - frequency Hz
            l - wavelength nm
            Å - wavelength Å
        norm : boolean, optional (default: False)
            normalize the lines to the total as determined by units
        split : boolean, optional (default: False)
            split into two plots based on atomic and molecular lines
        write : boolean, optional (default: False)
            write intensity matrix to logs 
        xlim : len-2 tuple of floats, optional (default: None)
            abscissae limits in cm**3 s**-1, defaults to default limits
        linewidth : float, optional (default: 3)
            width of lines plotted 

        Returns
        -------
        matplotlib.pyploy.Figure object
        '''
        from matplotlib.pyplot import figure

        try:
            ax = fig.get_axes()[0]
        except:
            fig = figure(figsize=figsize)
        if units == 'ev':
            xunit = r'E [eV]'
            if xlim is None:
                xlim = (0,20)
        elif units == 'v':
            xunit = r'v [$\mathrm{cm^{-1}}}$]'
            if xlim is None:
                xlim = (1.4e4,1.05e5)
        elif units == 'f':
            xunit = 'Frequency [THz]'
            if xlim is None:
                xlim = (400,3500)
        elif units == 'l':
            xunit = r'$\rm{\lambda}$ [nm]'
            if xlim is None:
                xlim = (80,700)
        elif units == 'Å':
            xunit = r'$\rm{\lambda}$ [Å]'
            if xlim is None:
                xlim = (800, 7000)
        if norm is True:
            yunit = 'Normalized intensity []'
        else:
            yunit = r'Counts [$\rm{s^{-1}}]$'



        if intensity is None:
            data = self.intensity(Te, ne, E=E, Ti=Ti, ni=ni, Srec=Srec, n0=n0,
                            units=units, norm=norm, write=False, **kwargs)
        else:
            data = intensity
        if split is True:
            species = ['atomic','molecular']
            for i in range(2):
                ax = fig.add_subplot(3,1,i+1)
                for d in range(len(data[i][0,:])):
                    if data[i][0, d] != 0:
                        ax.semilogy([data[i][0, d], data[i][0, d]],
                                [0, data[i][1, d]],'k', linewidth=linewidth)
                ax.set_ylim(0, 1.1*max(data[i][1,:]))
                ax.set_xlim(xlim)
                ax.set_xlabel(xunit)
                ax.set_ylabel(yunit)
                ax.set_title('Total radiation from {} processes'
                        ''.format(species[i]))
            ax = fig.add_subplot(313)
        else:
            ax = fig.add_subplot(111)
        color = ['b', 'r']
        if norm is True:
            data = self.intensity(Te, ne, Ti=None, ni=None, E=E, units=units,
                    norm=False, write=False)
            n = sum(data[0][1,:]) + sum(data[1][1,:])
        else:
            n = 1
        for i in range(2):
            ax.semilogy([], [], 'b-')
            ax.semilogy([], [], 'r-')
            for d in range(len(data[i][0,:])):
                if data[i][0, d] != 0:
                    ax.semilogy([data[i][0, d], data[i][0, d]], 
                            [0,data[i][1, d]/n], linewidth=linewidth,
                            color=color[i])
            try:
                ax.set_ylim((None, 1.1*max(max(data[0][1,:]), max(data[1][1,:]))/n))
            except: 
                pass
                #ax.set_ylim((0,None))
            ax.set_xlim(xlim)
            ax.set_xlabel(xunit)
            ax.set_ylabel(yunit)
            ax.set_title('Total radiation from molecular processes')
            #ax.legend(['Atomic bands', 'Molceular bands'])
        fig.show()


    def plotrate(
            self, database, h123, name, T, n, E=0.1, logx=True, logy=True, res=200, 
            color='k', linestyle='-', figsize=(12,13/1.618033), title='',
            ylim=[1e-14,1e-6], linewidth=2, savename=None, figtype='png',
            xlim=None, ncol=3, ax=None, **kwargs):
        ''' Plots the rate of the requested reaction

        Parameters
        ----------
        database : string
            the database where to look for the rate
        name : string
            the handle of the reaction
        T : len-2 tuple of floats
            the temperature boundaries, in eV, for which to plot the 
            rate. Assumed to be Te or Ti, depending on the reaction
        n : float 
            the densities, in cm**-3, for which to plot the rate. 
            Assumes to be ni or ne, depending on the reaction
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        logx : boolean, optional (default: True)
            plots a logarithmic x-axis
        logy : boolean, optional (default: True)
            plots a logarithmic y-axis
        res : int, optional (default: 200)
            the number of points to be sampled in the temperature 
            interval
        color : string, optional (default: 'k')
            the plot color
        linestyle : string, optional (default: '-')
            the plot line style
        fig : matplotlib.pyplot.figure, optional (default: None)
            figure which to plot on. By default creates a new figure and
            returns the figure handle
        figsize : tuple (float, float), optional (default: (14, 8.33))
            size of the created figure
        savename : string, optional (default: None)
            Name to save figure by. Figures are saved to 
            {Crumpet.path}/output{savename}.{figtype}. By default, 
            the figure is not saved.
        ylim : len-2 tuple of floats, optional (default: (1e-7, 1e7))
            ordinate limits in cm**3 s**-1
        linewidth : float, optional (default: 3)
            width of lines plotted 
        figtype : string, optional (default: 'png')
            filetype to save figure as. See matplotlib.pyplot.savefig
            for options
        title : string, optional (default: '')
            the figure title
        ncol : int, optional (default: 3)
            ncol keyword of legend: how many columns in legend
        xlim : len-2 tuple of floats, optional (default: None)
            abscissae limits in cm**3 s**-1, defaults to default limits

        Returns
        -------
        matplotlib.pyploy.Figure object
        '''
        from numpy import linspace, logspace, log10
        from matplotlib.pyplot import figure

        if ax is None:
            fig = figure(figsize=figsize)
            ax = fig.add_subplot(111)
        if logx is True:
            pf = ax.semilogx
            try:
                x = logspace(log10(T[0]), log10(T[1]),res)
                y = [self.get_rate(database, h123, name, Te=i,Ti=i, ne=n , ni=n, E=E, **kwargs) for i in x]
                xlabel = 'Temperature [eV]'.format(T)
                titapp = (', '*(len(title) > 0) + r'n={:.2E} '.format(n)
                        + r'$\rm{{cm^{{-3}}}}$')
            except:
                x = logspace(log10(T[0]), log10(T[1]),res)
                y = [self.get_rate(database, h123, name, Te=i, Ti=i, ne=n, ni=n, E=E, **kwargs) for i in x]
                xlabel = r'Density [$\rm{cm^{-3}}$]'
                titapp = ', '*(title != '') + r'T={} eV'.format(T)
            if logy is True:
                pf = ax.loglog
        else:
            pf = ax.plot
            try:
                x = linspace(T[0], T[1], res)
                y = [self.get_rate(database, h123, name, i, n, E) for i in x]
                xlabel = 'Temperature [eV]'.format(T)
                titapp= ', '*(len(title) > 0) + r'n={:.2E} '.format(n)+\
                        r'$\rm{{cm^{{-3}}}}$'
            except:
                x = linspace(T[0], T[1], res)
                y = [self.get_rate(database, h123, name, Te=i,Ti=i, ne=n , ni=n, E=E) for i in x]
                xlabel = r'Density [$\rm{cm^{-3}}$]'
                titapp = ', '*(title != '') + r'T={} eV'.format(T)
            if logy is True:
                pf = ax.semilogy
        pf(x, y, color=color, linestyle=linestyle, linewidth=linewidth,
                label=database + '.' + name)#, marker='.')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$\rm{Rate}$ [$\rm{cm^{3}/s}$]')
        ax.set_title(title + titapp)
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
            fig.savefig('output/figs/{}.{}'.format(savename, figtype), dpi=300,
                    edgecolor=None, format=figtype, bbox_inches='tight')
        try:
            fig.show()
            return fig
        except:
            pass
            return

    def write_reaction_index(self, outname='reaction_index_{}.csv', outpath='output'):
        from pandas import DataFrame

        df = DataFrame(columns=['Identifier', 'Database', 'Reaction type',
                        'Reactants', 'Products', 'Source file', 
                        'Scaling factor'])

        i = 0
        for database, dbdata in self.reactions.items():
            for rtype, rdata in dbdata.items():
                for reaction, readata in rdata.items():
                    if rtype.upper() == 'RELAXATION':
                        src = '{:.4E}'.format(readata.coeffs)
                    else:
                        try:
                            src = readata.sourcefile.split('/')
                            if ('H2' in src[-2]) or ('D2' in src[-2]):                        
                                src = '{}/{}'.format(src[-2], src[-1])
                            else:
                                src = '{}/{}/{}'.format(src[-3], src[-2], src[-1])
                        except:
                            src = readata.sourcefile.split('/')[-1]

                    buff = readata.print_reaction().split(':')[1].strip()
                    df.loc[i] = [ reaction, database, rtype,
                            buff.split('=>')[0].strip(),
                            buff.split('=>')[1].strip(),
                            src, 
                            ''*(readata.scale==1)+str(readata.scale)*(readata.scale!=1)]
                    i += 1
        try:
            mkdir('{}/{}'.format(self.path, outpath))
        except:
            pass
        df.sort_values('Identifier')
        df.to_csv('{}/{}/{}'.format(self.path, outpath, 
                outname.format(self.isotope)))


    def lifetimes(self, Te, ne, E=0.1, Ti=None, ni=None, Srec=True, n=None):
        ''' **INCOMPLETE** Species lifetimes perturbation from equilibrium

        Retrurns a ordered list (according to the input file species) of 
        the species lifetimes. The species lifetimes are calculated by 
        an intial guess and a set of sink/source terms that are evolved
        to a steady-state. Each species of the steady-state distribution
        is then perturbed in order, and a exponential fit is applied to 
        the decay of the perturbation, determining the lifetime.

        Parameters
        ----------
        Te : float
            electron temperature in J
        ne : float
            electron density in m**-3
        Ti : float, optional (default: None)
            ion temperature in eV. Default sets Ti=Te
        ni : float, optional (default: None)
            ion density in cm**-3. Default sets ni=ne
        E : float, optional (default: E=0.1)
            molecular energy in eV for rates
        Srec : boolean, optional (default: True)
            switch for considering 'external' sinks/sources (e.g. 
            re-association)
        n : list of floats, optional (default: None)
            an initial guess for the steady-state densities of length N. 
            Defaults to n0 defined in the input

        Returns
        -------
        list of floats
            Radiative lifetimes of each species as calculated by fitting
            an exponential to the decay of the perturbation
        '''
        from numpy import log, ones, matmul, polyfit, zeros, exp, sum
        from scipy.optimize import fsolve, minimize
        from scipy.linalg import norm
        from scipy.integrate import solve_ivp
        from matplotlib.pyplot import figure

        def dt(t, n, mat, ext, n0):
            return matmul(mat, n) + ext

        def findfac(fac, n, mat, ext, n0):
            nend_sum = self.totpart(solve_ivp(lambda x, y: dt(
                    x, y, mat, ext, n, fac), (0, 1), n, 'LSODA',
                    dense_output=True).y[:,-1])
            return abs(self.totpart(n0) - nend_sum)

        na = 3.23058E+15
        nm = 1.81807E+16
        diva = -1.10649E+19
        divm = -2.45590E+11
        srca = 0.00000E+00
        srcm = 6.24142E+18
        bga = 8.34627E+09
        bgm = 4.63451E+10
        psor = -1.99728E+17
        psorgc = 6.25988E+14
        vol = 1e4
        # TODO: Check that volumes go right: CM to CM
        # TODO: Implement approximation of pumping 
        # Use input n0 as default unless explict n0 requested
        if n is None: 
            n=self.n0() 
        #NN=len(self.species)*len(self.species)
        n[0] = 1e-6*na
        n[1] = 1e-6*nm
        # Get the full rate matrix
        mat, ext = self.M(Te, ne, Ti, ni, E, write=False) 
        # TODO: use Greenland CRM to assess the SS - the same assumptions
        ''' Assume re-association to achieve SS '''
        # TODO: fix ext sources and dissociation
        print('REVISE')
        #mat[0,0]=-reassociation
        #mat[1,0]=-mat[0,0]*0.5 # Conserve nucleii
        ext[0] = (psorgc + diva + srca + bga)
        ext[1] = (divm + srcm + bgm)
        mat[0,0] = (psor/ne)
        # Get matrices and eigenvalues/-vectors
        

        M, T, D = self.gl_crm(mat, ext, matrices=True, Srec=Srec) 
        # Get matrices and eigenvalues/-vectors
        Meff, extp, n0p = self.gl_crm(mat, ext, n=n) 
        for i in D:
            print('{:.3E}'.format(-1/i))
        ss_sol = solve_ivp(lambda x, y: dt(x, y, M, ext, n),
                (0,10), n, 'LSODA', dense_output=True)
        f = figure(figsize=(7, 7))
        ax = f.add_subplot(111)
        for i in range(len(n)):
            line = '-' + '-'*('H2' in self.slist[i])
            ax.semilogx(ss_sol.t, ss_sol.y[i,:], line, label=self.slist[i]) 
            if i == 0:
                ax.semilogx(ss_sol.t, ss_sol.y[i,:], 'k.', label=self.slist[i])
        ax.legend(ncol=3)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Density [cm**-3]')
        print('UEDGE molecular content={:.5E} m**-3'.format(nm))
        print('CRUMPET molecular content={:.5E} '
                'm**-3'.format(1e6*self.totmol(ss_sol.y[:,-1])))
        return
        
        global ss_n
        ss_n = ss_sol.y[:,-1]
        ret = [1e20]
        # Perturb each species
        for i in range(1, len(n)):
            try:
                pert = zeros((len(n),))
                pert[i] = 1e6
                sol = solve_ivp(lambda x, y: dt(x, y, mat, ext, n), (0, 1e-3),
                        ss_n + pert, 'LSODA', dense_output=True, events=ss)
                a, b = polyfit(sol.t, log(sol.y[i,]), 1)
                print('tau({})'.format(self.slist[i]).ljust(20,' ')+
                        '= {:.3E}'.format(-1/a))
            except: 
                pert = zeros((len(n),))
                pert[i] = 1e6
                sol = solve_ivp(lambda x, y: dt(x, y, mat, ext, n), (0, 1e-3),
                        ss_n + pert, 'LSODA', dense_output=True, events=ss)
                a, b = polyfit(sol.t, log(sol.y[i,]), 1)
                print('tau({})'.format(self.slist[i]).ljust(20,' ')+
                        '= {:.3E}'.format(-1/a))
            ret.append(-1/a)
                #sol=solve_ivp(lambda x,y: dt(x,y,mat,ext),(0,1e-3),ss_n+pert,
                #                'LSODA',dense_output=True)
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
        for i in range(1, len(n)):
            print('tau({})'.format(self.slist[i]).ljust(20,' ')+
                    '= {:.3E}'.format(-1/mat[i,i]))
        return ret 


