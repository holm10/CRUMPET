=========
Reactions
=========

The reaction card defines reactions included in the CRM that are read from one of the default databases defined under the **rates** card. The subcards define the database where to look for the reaction, and the reaction name/number/id, separated by an underscore: :code:`database_ID`. Here, :code:`database` is one of the databases defined in **rates**, and :code:`ID` is one of the reactions available in the database. The following built-in databases are also available to the user, without setting up them in **rates**: for these database, ID is only an identifying handle.

- JOHNSON: Einstein A-coefficients between electronic levels in hydrogen atoms. ioffe.ru/astro/QC/CMBR/dscrpt.pdf
- APID: Ionization rates from electronically excited atom levels. R. K. Janev and J. J. Smith, Cross Sections for Collisions Processes of Hydrogen Atoms with Electrons, Protons and Multiply Charged Ions, vol. 4. Vienna: International Atomic Energy Agency, 1993.

The Following line **must** define the reaction using the handles defined in **species** and **background**. The syntax is :code:`R1 + R2 > X*P1 + Y*P2 + ...`. Here, :code:` > ` separates the reactants (LHS) from the products/fragments (RHS). :code:`R1` and :code:`R2` are the reactant labels, one of which belongs to **background**, and one of which belongs to **species**: a maximum of two reactants can be specified in the present implementation, but only one reactant is possible (decay, etc). :code:`P1`, :code:`P2`, ... defines the products: there are no limitations on the reaction products. Here, the multipliers :code:`X` and :code:`Y` are used as shorthand notations: they do not have to be integers, in order to allow for statistical distributions. However, it is up to the user to supply reactions that conserve particles. Note that only collisional reactions with one background species and one CRM species on the LHS can be defined.

Each reaction subcard can also be assigned the mean kinetic energy of the products on the following lines. By default, the mean kinetic energy of the products is 0. The syntax for the mean energy of the products is :code:`K=X`, where "X" is either a float, but can also be an expression of simple binary operators, e.g. :code:`2*3+4`. The mean kinetic energy is extracted from the background reactant (electron or ion/atom) and given to the products: in case of proton-impact reactions, the kinetic energy is given to the electrons. The mean kinetic energy can also contain the strings "Te", "Ti", or "Tm": this substitutes the local electron, ion, or molecular temperature, respectively, as the kinetic energy of the products. This way, the electrons in an ionization reaction can be defined to be heated to the local electron temperature by the reactant.

The special characters :code:`$` and :code:`&` can be utilized in excitation/de-excitation reactions for the electronic and vibrational levels. Here, the electronic or vibrational levels, up to :code:`nmax` or :code:`vmax` as defined in **settings** if the strings "n=" or "v=" are included in the species labels, respectively. The same substiution is made in :code:`ID`, so the same symbols must be included in :code:`ID` in order to create electronically/vibrationally dependent rates. Here, the hierarchy between :code:`$` and :code:`&` is such, that :code:`$` is always evaluated first, and :code:`&` is only evaluated if :code:`$` is called. Note that, presently, both :code:`$` and :code:`&` must be applied to the same processe, i.e. either electronic or vibrational processee, but not separately to each.

The special subcard :code:`* CUSTOM` can be used to define a set of user-specified rates, as explained in **custom**. Instead of a reaction, :code:`CUSTOM` reactions must be followed by the relative path to the custom rate data file, relative to Crm.path.

An example definition of the REACTION card:

.. code-block::

	** REACTIONS
	* HYDHEL_2.2.12
	e + H2+ > e + p + H(n=1)
		K=2*4.3
	* H2VIBR_2.$v&
	e + H2(v=$) > e + H2(v=&)
	* CUSTOM
	rates/custom.dat

.. toctree::

.. toctree::
    :hidden:
    
    custom
