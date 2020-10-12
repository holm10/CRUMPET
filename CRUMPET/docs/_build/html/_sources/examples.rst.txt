CRUMPET examples
================

Example input file
******************

.. include:: input.txt


Example custom reactions
************************

.. include:: custom.txt


Generate a CRM form the input file
**********************************

.. code-block::
	
	import CRUMPET
	crm = CRUMPET.Crumpet('input.dat') # Create a CRUMPET object, defined in input.dat

Time-dependently evolve the initial distribution
************************************************


.. code-block::
	
	t_stop = 1e-4 # Stop-time for the time-evolution in s
	Te = 5 # Electron temperature in eV
	Ti = 5.5 # Ion temperature in eV
	ne = 1e13 # Electron density in cm**-3
	ni = 2e13 # Ion density in cm**-3
	conservation = True # Show particle and energy conservation
	gl = False # Don't evaluate the Greenland CRM: solve the full set of ODEs
	Qres = True # Plot the Q-species

	crm.plot_crm(t, Te, ne, Ti=Ti, ni=ni, gl=gl, conservation=conservation, Qres=Qres)

Plot a synthetic spectra
************************


.. code-block::
	
	Te = 5 # Electron temperature in eV
	Ti = 5.5 # Ion temperature in eV
	ne = 1e13 # Electron density in cm**-3
	ni = 2e13 # Ion density in cm**-3

	from numpy import zeros
	ext = zeros((len(crm.species),)) # Initialize an array for the external sink/sources of each species in the CRM
	ext[:2] = [-2e11, 1e12] # Define an external atom sink and external molecular source (setup-specific species)
	ionization = -1e8 # Ionization sink of atoms (setup-specific species)

	crm.spectrum(Te, ne, Ti=Ti, ni=ni, ext=ext, ionization=ionization)


