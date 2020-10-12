======================
The CRUMPET input file
======================

Crumpet uses an input file using a custom syntax, outlined here. The CRUMPET input file consists of a number of "cards", setting up the main parts of the code, followed by a number of "subcards", setting up the cards. The cards are the top-level entries, and are marked by ``** ``, whereas the cards are marked by ``* ``, both followed by a space and the appropriate entry. The setup of each subcard is explained in detail below. The following are universally true for CRUMPET:

- Case insensitive syntax
- `#` mark comments
- Indentation levels are not used
- The input file is read line-by-line


The CRUMPET input file consists of a number of subards, defining the CRM setup:

.. toctree::
 
.. toctree::
    
    species
    background
    reactions
    rates
    settings

