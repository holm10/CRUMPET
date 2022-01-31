=======
Species
=======

The species card contains subcards for each species to be included in the CRM model. Each species must be given a **unique** handle, used for identifying the species in the code execution. The only options for the subcards is the potential level of the species: the potential energy level is only relative to the other species and the background. The potential is defined on the line below the subcard and handle as "V X", where "V" identifies the setting, and X is attributed as the species potential. Note that the only accepted separators between the label and the value are whitespeces, at least one of which are required.

The species label **cannot** contain any whitespaces, and the molecular species are identified as any label containing the strin "{isotope}2", where {isotope} is defined in the **settings** subcard, and the atom species are identified as any species containing the string "{isotope}" without the "2". The strings "n=X" and "v=X", where "X" is an integer are reserved for marking the electronic levels of atoms and the vibrational levels of molecules, respectively. Any species label including the string "n=X" is treated as an atom at electron level "X". Analogously, any species label including the string "v=X" is treated as a molecule at vibrational level "X". The strings "e" and "p" are reserved for background species.

An example of defining ground-state atoms and molecules, along with ionic ions are defined below. Here, the zero-potential is defined as the minimum potential of the ground state as a function of intermolecular separation. 

.. code-block::

	# Define the species to be included
	** SPECIES # Begin block
	* H(n=1) # Ground-state atoms
		V 2.375
	* H2(v=0) # Ground-state molecules
		V 0.27504
	* H2+ # Ionic molecules
		V 15.56

Note that the CRM is only valid for collisional reactions involving one CRM species and one background species.

.. toctree::

.. toctree::
    :hidden:
    :maxdepth: 3

