Collisional-Radiative theory
============================
Here, a short overview of collisional-radiative (CR) theory is given to provide context for the CRUMPET code. The documentation is based on [1] and [2], where a complete description of the CR model implemented can be found. 


The Greenland CRM equation
**************************
The time-evolution of neutral species interacting with a background plasma is described by a set of nonlinear rate equations of the form 

.. math::
    
    \dot{n_k}=\sum_{i,j}\mathcal{R}_{i,j}^k(T)n_in_j+\sum_j A_j^kn_j-n_k\sum_{i,j}\mathcal{R}_{i,k}^j(T)n_i-n_k\sum_jA_k^j+\Gamma_k,

where :math:`n_k` is the density of species :math:`k`, :math:`\mathcal{R}_{i,j}^k(T)` is the rate coefficient for production of species :math:`k` from collisions between species :math:`i` and :math:`j` at plasma temperature :math:`T,~A_j^k` is the Einstein A coefficient for radiative decay from :math:`j\rightarrow k`, and :math:`\Gamma_k` is an external sink/source of species :math:`k`, such as recombination or transport.

The set of equations can be linearized by fixing appropriate densities :math:`n_1\cdots n_m`, specifying the background plasma environment:

.. math::

    \dot{\mathbf{n}}=\mathbf{M}(n_1\cdots n_m)\mathbf{n}+\mathbf{\Gamma}(n_1\cdots n_m).

Here, :math:`\mathbf{M}(n_1\cdots n_m)` is the full :math:`N\times N` matrix of rates at plasma temperature :math:`T_e` and :math:`T_i` and :math:`\mathbf{\Gamma}(n_1\cdots n_m)` is a source/sink provided by collisions between pairs of backgrounds species, and any external sinks and sources; :math:`\bf n` is a vector of the states whose dynamics are evaluated.

The states to be evaluated are subsequently divided into two classes, *P* and *Q*, which contain the species to be evaluated time-dependently and the conduit species, respectively. The corresponding projection operators *P* and *Q* can then be defined: by ordering the states *k* in the linearized equations, the equation can be expressed as

.. math::

    \begin{bmatrix} \dot{\mathbf{n}}_P \\ \dot{\mathbf{n}}_Q \end{bmatrix} = \begin{bmatrix} \mathbf{M}_P & \bf{H} \\ \bf{V} & \mathbf{M}_Q \end{bmatrix} \begin{bmatrix} \mathbf{n}_P \\ \mathbf{n}_Q \end{bmatrix}	+ \begin{bmatrix} \mathbf{\Gamma}_P \\ \mathbf{P}_Q \end{bmatrix},

where :math:`P\mathbf{M}P=\mathbf{M}_P`, :math:`P\mathbf{M}Q=\mathbf{H}`, :math:`Q\mathbf{M}P=\mathbf{V}`, :math:`Q\mathbf{M}Q=\mathbf{M}_Q`, :math:`P\mathbf{n}=\mathbf{n}_P`,  :math:`Q\mathbf{n}=\mathbf{n}_Q`, :math:`P\mathbf{\Gamma}=\mathbf{\Gamma}_P`, and :math:`Q\mathbf{\Gamma}=\mathbf{\Gamma}_Q`. Here, there are :math:`N_P` *P*-states and :math:`N_Q` *Q*-states.

It has been shown that the solution of the P-space part of the full, linearized set of equations can be reproduced by the following equation: 

.. math::

 	\dot{\mathbf{n}}=\mathbf{M}_{eff}\mathbf{n}'_P+\mathbf{\Gamma}'_P.

The validity of the above equation is determined using the eigenvectors of **M**: we define **T** as the matrix whose columns are the normalized eigenvectors of **M** and :math:`\mathbf{D}=\mathrm{diag}\{\lambda^{(1)},\lambda^{(2)},\ldots,\lambda^{(N)}\}` in order of increasing order of the absolute magnitude of the eigenvalues :math:`\lambda` such that

.. math::
	
	\bf MT = TD

holds. Analogous to above, the projection operators can be used to subdivide **T**:

.. math::

	\begin{bmatrix} \mathbf{T}_P & \bf{\Delta} \\ \bf{\delta} & \mathbf{P}_Q \end{bmatrix}.

The validity of the CRM equation, thus, requires :math:`\delta` and :math:`T_Q^{-1}\delta\ll1`, where the equivalence refers to the largest Manhattan norm along the columns of the matrix. 

The remaining terms of the CRM equation are, hence, given by

.. math::

 	\mathbf{M}_{eff} = \mathbf{M}_P - \mathbf{HM}_Q^{-1}\mathbf{V}

and 

.. math::

	\mathbf{\Gamma}'_P = \mathbf{\Gamma}_P - \mathbf{\Delta T}_Q^{-1}\mathbf{\Gamma}_Q,

with initial conditions given by

.. math::
	
	\mathbf{n}'_P(0)=\mathbf{n}_P(0)-\Delta T_Q^{-1}\mathbf{n}_Q(0).


Model limitations/assumptions
-----------------------------
This CRM desctiption of the denity rate coefficients has the following limitations, as described in [1] and [2]:

- Only collisional reactions between a background species and a CRM species are considered: reactions between the CR-modeled species are not treated by this formulation


The CR energy rates
*******************
The energy rate coefficients :math:`S_j^k(T)` associated with a CRM can be constructed in analogy with the CRM equation above, as explained in [2] for the electrons. CRUMPET makes an attempt to extend the energy rate coefficients to other species/processes than the electrons, in a manner analogous to that outlined in [2]. In CRUMPET, the potential energy sink per unit volume, :math:`\mathcal{E}`, is solved from the equation

.. math::
	
	\dot{\mathcal{E}}_V^k=\sum_{i,j}\left( V_p - V_r \right)\mathcal{R}_{i,j}^k(T)n_in_j - n_k\sum_j(V_k-V_j)A_k^j,

where :math:`V_r` and :math:`V_p` refers to the total potential of the reactants and products, respectively, and :math:`V_k` and :math:`V_j` to the potential of CRM species *k* and *j*, respectively, and *i* is a background (non-CRM) species. The energy lost as radiation is solved from

.. math::
	
	\dot{\mathcal{E}}_{\gamma}^k=n_k\sum_j\left( V_k - V_j\right),

which is further divided into atom and molecule radiation depending on whether *j* and *k* are atomic or molecular species. For the Franck-Condon energy, CRUMPET applies the UEDGE formalism of a common energy and, thus, temperature for the ions and atoms. Hence, a common Franck-Condon sink/source term is used for the ion/atom energy:

.. math::
	
		\dot{\mathcal{E}}_{FC,i/a}^k=\sum_{i,j}\Delta_{i,j}^k\mathcal{R}_{i,j}^k(T)n_in_j,

where :math:`\Delta_{i,j}^k` is the mean kinetic energy of the reaction products, *K*, for reactions where electrons are one of the reactants, such as dissociation energy. 
.. If the reaction is a proton-impact reaction, i.e. one of the reactions is a plasma ion, producing an electron, :math:`\Delta_{i,j}^k=-T_e` is assumed, where :math:`T_e` is the local electron temperature. 
Here, the Franck-Condon energy resulting from production of conduit species is directly attributed to the ion/atom equation: thus, the energy of the conduit species are not evolved. Subsequently, any energy-dependent reaction rates (heavy-particle interactions) of reactions involving the conduit species are evaluated at the prescribed target energy (or, alternatively, temperature). This approximation is valid unless there is an energy source for the conduit species capable of significantly affecting the target particle energy-dependent reaction rates. Finally, :math:`\dot{\mathcal{E}}_e` is calculated as the sum of the above energy sinks/sources by assuming energy conservation of the system. 

The CRM evaluation is then done in analogue to the density rate coefficients:

.. math::

	\dot{\mathbf{\mathcal{E}}}=\mathbf{U}(n_1\cdots n_m)\mathbf{n}+\mathbf{\Xi}_{E}(n_1\cdots n_m),

where 

.. math:: 

	U_{k,j}=S_j^k(T)ne

and 

.. math::
	
	\Xi_{E,k}=\sum_i S_i^kn_en_i,

where :math:`S_j^k` is the appropriate term according to the above equations, *j* and *k* are CRM species, and *i* are background species. Analogous to above, the energy rate matrix **U**, ordered by *P* and *Q* space, can be expressed as a block matrix, using the previously defined projection operators:

.. math::
	
	\mathbf{U}=\begin{bmatrix} \mathbf{U}_P & \mathbf{U}_H \\ \mathbf{U}_V & \mathbf{U}_Q \\ \end{bmatrix}.

The subsequent CR approximation for the energy terms, thus, becomes:

.. math::
	
	\mathbf{\dot{\mathcal{E}}}_P=\mathbf{U}_{eff}\mathbf{n+\Xi}_E,

where 

.. math::

	\mathbf{U}_eff=\begin{bmatrix} 	\mathbf{U}_P - \mathbf{U}_H\mathbf{M}_Q^{-1}\mathbf{V} \\ %
									\mathbf{U}_V - \mathbf{U}_Q\mathbf{M}_Q^{-1}\mathbf{V} \\ \end{bmatrix},

where :math:`\mathbf{M}_Q` and :math:`\mathbf{V}` are the same as defined above. This yields an :math:`N\times N_P` matrix, rather than than the :math:`N_P\times N_P` matrix for the density rate coeffients: this reflects the fact that although a small population change is associated with the conduit species, a large energy change can still be associated with a conduit species if :math:`S_j^k` is sufficiently large. The above CRM is evaluated separately for each energy sink/source, i.e. for the electron energy sink/source, the ion/atom energy sink/source, the potential energy sink/source, the atom line radiation sink/source, and the molecule line radiation sink/source.

For coupling of CRUMPET to external codes, :math:`\mathbf{\dot{\mathcal{E}}}_P` is summed over the columns to yield the net energy change associated with the *P*-species densities, as the coupled code is assumed to only solve the densities of the *P*-species.


Model limitations/assumptions
-----------------------------
This CRUMPET description of the energy sinks and sources conforms to the UEDGE formalism and, presently, has the following limitations:

- Presently, the model evaluates a common ion/atom source, which is compatible with the UEDGE formalism
- The Franck-Condon energy is assumed to go directly into the ion/atom energy, rather than into any resulting conduit species: this assumption is valid unless the energy source is sufficiently strong to significantly alter the target particle energy-dependent reaction rates
..- Electrons stemming from proton-impact reactions are assumed to be heated to the local electron temperature by the impacting particle
- Electron-impact ionization reactions only assess the energy change associated with the ionization potential: any volumetric cooling due to the newly formed ion/electron pair is assumed to be considered by the coupled transport code
- The energy sink/source due to electron-ion recombination (radiative and three-body recombination) have not yet been implemented
- The energy sink/source from molecular re-association has not yet been implemented: considerations of the formulation of the CR formulation raise question about how to consider re-association (see above)


References
**********
| [1] P. T. Greenland, “Collisional–radiative models with molecules,” Proc. R. Soc. Lond. A, vol. 457, pp. 1821–1839, 2001.
| [2] P. T. Greenland, “The CRMOL Manual: Collisional Radiative Models for Molecular Hydrogen in Plasmas”. Forschnungszentrum Jülich report, Jül-3858, 2001.
