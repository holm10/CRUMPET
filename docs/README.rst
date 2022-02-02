CRUMPET: Collisional-Radiative UEDGE Model for Plasma-Edge Theory
======================================================================

Full manual at: https://holm10.github.io/CRUMPET/

What is CRUMPET?
----------------
CRUMPET is a versatile Python script for creating collisional-radiative models, originally developed for fusion plasma application, with the multi-fluid plasma-edge code UEDGE [1]. Collisional-radiative models are necessary under conditions where a system is under neither Saha/LT equilibrium, nor Coronal equilibrium. Under such conditions, the collisional and radiative processes of the system are competing, and the net effect of all such processes, which are often competing, need to be carefully evaluated. CRUMPET has the capability to consider molecular-assisted processes, such as recombination, dissociation, and ionization (MAR, MAD, and MAI, respectively) as well as Franck-Condon and potential energy sinks/sources and molecular and atomic line radiation.

Neutral dynamics involve a large set of electronically, vibrationally, and rotationally excited species, and a number of reactions associated with each such species. This creates a large number of species, and transitions if each species are to be tracked time-dependently. However, a large number of the species are so-called "conduit" species: their lifetimes are too short for any significant population to survive in the plasma. Rather, such conduit species mediate further reactions yielding more stable particles, which exist in significant populations in the plasma. Simulating all species and reactions, including conduit-species, time-dependently requires vast computational resources. Instead, a collisional-radiative model can be applied to the system, yielding effective reaction rates for only a subset of species which are to be simulated time-dependently. 

Atomic CRMs have been widely applied to fusion plasmas for the last decades as they are straightforward to formulate. However, molecular CRMs impose additional conditions on the formulation of the CRM compared to atomic CRMs. These conditions are outlined in [2], and have been implemented in CRUMPET in order to create molecular collisional-radiative models.

Features
--------
- Versatile creation of collisonal-radiative models
- CRM formulations only dependent on data supplied to the model
- Not bound to a specific element or isotope
- Capability to read rate data from ADAS and EIRENE data files
- Capability to use custom rate data
- Solve the full system of ODEs
- Create a Greenland-style CRM [2] without any underlying assumptions of the background plasma
- Solve the energy sinks and sources for a number of processes
- Write rate coefficients to tabulated data files
- Assess the accuracy of the current CRM
- Use Greenland's algorithm [2] to find the optimal CRM for a given plasmas temperature and density
- Solve the steady-state CRM populations using external sinks and sources (e.g. from transport codes)
- Create synthetic spectra for steady-state plasmas, which can be used for post-processing plasma-edge simulations 

How is it used?
---------------
Crumpet reads a CRUMPET-specific input file where the CRM is set up. The CRM species, background species, and time-dependent/conduit species are defined in the input file, as are the reactions to be considered. The reactions use rate data supplied by the user: a nnumber of standard A&M databases, such as ADAS and EIRENE can also be read in. A number of options and settings can also be used to control the behavior and setup of the CRM. 

Installation
------------
CRUMPET is installed by downloading the git repository to either a location that is listed in PYTHONPATH, or by adding the path to the downloaded repository to PYTHONPATH.

Support
-------
If you have issues, or would like to use CRUMPET for a certain task and need assistance, you can reach the author via andreas.holm{at}aalto.fi.

Contributions
-------------
Anyone is free to contribute to CRUMPET via github using personal forks or branches with a merge-request.
- Issue tracker: github.com/holm10/CRUMPET/issues
- Source code: github.com/holm10/CRUMPET

License
-------
CRUMPET is licensed under the MIT License 

References
**********
| [1] T. D. Rognlien and M. E. Rensink, "Users manual for the UEDGE edge-plasma transport code", LLNL Report, 2017. https://github.com/LLNL/UEDGE
| [2] P. T. Greenland, “Collisional–radiative models with molecules,” Proc. R. Soc. Lond. A, vol. 457, pp. 1821–1839, 2001.
