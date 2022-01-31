========
Settings
========

The **settings** card sets up the CRM. **settings** can be given the following subcards:

- NP: The P-space dimension. The P-species are chosen as the NP first species defined in **species**
	- :code:`* NP X`, where X is an integer (default: 2)
- ISOTOPE: The character or string sequence identifying the CRM isotope. See **species** for further information
	- :code:`* ISOTOPE X`, where X is a string (default: H)
- MASS: The isotope mass in AMU. Presently only used for the custom rate SIGMA
	- :code:`* MASS X`, where X is a float (default: 1)
- VMAX: The maximum vibrational level of molecules to be considered. Considers vibrational levels in the interval [0,VMAX]
	- :code:`* VMAX X`, where X is an integer (default: 14)
- NMAX: The maximum electronic level of atoms to be cosidered. Considers electronic levels in the interval [1,NMAX]
	- :code:`* NMAX X`, where X is an integer (default: 8)
- VERBOSE: Switch whether to output the CRM creation diagnostic to the prompt. 
	- :code:`* VERBOSE X`, where X is 0/1 or True/False (default: False)
- N0: The initial distribution of the CRM species :code:`* ISOTOPE` followed by handles and densities
	- :code:`handle X`, where :code:`handle` is one of the handles defined in **species** and :code:`X` is a float of the density in cm**-3

An example of a SETTING card is given below:

.. code-block::
	** SETTINGS
	* vmax      14
	* nmax      8
	* Np        2
	* Isotope	H
	* Mass		1
	* n0
		H(n=1)  0e12
		H2(v=0) 1e10
	* verbose   0   


.. toctree::

.. toctree::
    :hidden:
