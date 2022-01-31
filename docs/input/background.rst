===========
Background
===========

The background species follow the same syntax as the CRM species. However, the background species must at least include electrons and plasma protons, marked by "e" and "p", respectively. The electrons should not have any associated potential, whereas the plasma proton potential must be calculated relative to the same level as in species. Further species beyond the plasma electrons and protons may be included. Note that the CRM is only valid for collisional reactions involving one CRM species and one background species.

An example background species card setup:

.. code-block::

	** BACKGROUND
	* e
	* p
		V 15.975


.. toctree::

.. toctree::
