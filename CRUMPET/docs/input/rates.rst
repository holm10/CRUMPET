=====
Rates
=====

Here, default databases can be defined to be read. **reactions** can use these databases to access predefined reaction rates. **rates** take no subcards, only parameters in the form :code:`database path`, where :code:`database` defines the database to be read, and :code:`path` the path to the source file relative to Crm.path. The following datbases are available to be read:

	- AMJUEL: The .tex file of the EIRENE A&M datbase AMJUEL
	- HYDHEL: The .tex file of the EIRENE A&M datbase HYDHEL
	- H2VIBR: The .tex file of the EIRENE A&M datbase H2VIBR
	- UE: The DEGAS2-rate file ehr*.dat
	- ADAS: ADAS ADF04 database

A RATE card example:

.. code-block::
	** RATES
		AMJUEL  rates/amjuel.tex # Inline comments now work too!
		HYDHEL  rates/hydhel.tex
		H2VIBR  rates/h2vibr.tex

.. toctree::

.. toctree::
    :hidden:

