=================
CUSTOM rate files
=================

Custom rate files use the a similar syntax as **reactions**, with the following changes:

- Only :code:`$` is evaluated
- Automated evaluation of :code:`$` only available for vibrational levels
- The card is name is used as the database handle 
- The subcards are of the form :code:`* TYPE ID`. Here, "TYPE" defines one of the following reaction rate types:
	- RATE: a one-dimensional polynomial fit, for the EIRENE data
		- 9-polynomial fit in rising order
	- COEFFICIENT: a decay coefficient in s**-1, similar to the Einstein A coefficient
		- A single float
	- SIGMA: a fit for the reaction cross-section, which is integrated over a Maxwellian to yield the cross-section. Fit taken from K. Sawada and T. Fujimoto, “Effective ionization and dissociation rate coefficients of molecular hydrogen in plasma,” Journal of Applied Physics, vol. 78, no. 5, pp. 2913–2924, 1995, doi: 10.1063/1.360037.
		- Eth, q0, A, Omega, W, gamma, nu: as in the paper appendix

Just as for **reactions**, the line immediately after the subcard **must** contain the reaction definition, using the same syntax. The following line must then contain the appropriate reaction coefficients for the specified reaction type. If :code:`$` is utilized, the line after the reaction definition must contain a "v", followed by the rate coefficients in the appropriate form on the following line. The "v"-coefficient pair must be repeated on the following lines :code:`vmax` in order for all the required rates to be included in the model. 

Custom rate file example: 

.. code-block::
	
	** CUSTOMRATES
	* COEFFICIENT c-diss
	H2(n=c) > 2*H(n=1)
		1E+05

	* COEFFICIENT B-$
	H2(n=B) > H2(v=$)

		# v 	coefficient
		v= 0		
		8.2772736E+07
		v= 1		
		7.5957320E+07
		v= 2		
		6.8873384E+07
		v= 3		
		6.7167080E+07
		v= 4		
		6.3768032E+07
		v= 5		
		6.3601840E+07
		v= 6		
		5.9270400E+07
		v= 7		
		5.7115708E+07
		v= 8		
		5.3714124E+07
		v= 9		
		5.1850780E+07
		v= 10		
		4.8581508E+07
		v= 11		
		4.6469164E+07
		v= 12		
		4.4187192E+07
		v= 13		
		3.8139332E+07
		v= 14		
		2.2083192E+07

	* RATE $-B
		e + H2(v=$) > e + H2(n=B)
			#	b0		b1		b2		b3		b4		b5		b6		b7		b8
			v= 0
			-0.2880946E+02	0.9541977E+01	-0.3966258E+01	0.1080775E+01	-0.2065085E+00	0.2690762E-01	-0.2253846E-02	0.1086306E-03	-0.2279712E-05
			v= 1
			-0.2816380E+02	0.9400263E+01	-0.4170421E+01	0.1253500E+01	-0.2667564E+00	0.3837598E-01	-0.3494823E-02	0.1803441E-03	-0.3999552E-05
			v= 2
			-0.2646792E+02	0.1037208E+02	-0.6565566E+01	0.2694575E+01	-0.7057833E+00	0.1153278E+00	-0.1135514E-01	0.6156517E-03	-0.1410782E-04
			v= 3
			-0.2700517E+02	0.8874048E+01	-0.4168608E+01	0.1355566E+01	-0.3117252E+00	0.4779701E-01	-0.4572125E-02	0.2449583E-03	-0.5592066E-05
			v= 4
			-0.2689817E+02	0.9289897E+01	-0.4689720E+01	0.1663134E+01	-0.4132763E+00	0.6734438E-01	-0.6745774E-02	0.3743037E-03	-0.8779752E-05
			v= 5
			-0.2599713E+02	0.8826899E+01	-0.4892533E+01	0.1934503E+01	-0.5211576E+00	0.8940473E-01	-0.9249642E-02	0.5241729E-03	-0.1247418E-04
			v= 6
			-0.2584816E+02	0.8795786E+01	-0.4832247E+01	0.1885146E+01	-0.5028363E+00	0.8574919E-01	-0.8841014E-02	0.5000067E-03	-0.1188453E-04
			v= 7
			-0.2570931E+02	0.8763256E+01	-0.4780183E+01	0.1843632E+01	-0.4873893E+00	0.8264650E-01	-0.8491882E-02	0.4792400E-03	-0.1137539E-04
			v= 8
			-0.2562608E+02	0.8767379E+01	-0.4755383E+01	0.1816035E+01	-0.4762224E+00	0.8032676E-01	-0.8226280E-02	0.4632726E-03	-0.1098098E-04
			v= 9
			-0.2554579E+02	0.8742564E+01	-0.4722188E+01	0.1790030E+01	-0.4665302E+00	0.7837564E-01	-0.8006431E-02	0.4501871E-03	-0.1066007E-04
			v= 10
			-0.2470037E+02	0.7878417E+01	-0.4423302E+01	0.1787829E+01	-0.4930125E+00	0.8625272E-01	-0.9064225E-02	0.5201051E-03	-0.1250256E-04
			v= 11
			-0.2467941E+02	0.7636692E+01	-0.4073292E+01	0.1562898E+01	-0.4164065E+00	0.7138478E-01	-0.7411070E-02	0.4220068E-03	-0.1009312E-04
			v= 12
			-0.2471478E+02	0.7599960E+01	-0.4091470E+01	0.1588046E+01	-0.4265000E+00	0.7344852E-01	-0.7644017E-02	0.4358356E-03	-0.1043086E-04
			v= 13
			-0.2490949E+02	0.7596427E+01	-0.4140032E+01	0.1632552E+01	-0.4431187E+00	0.7674929E-01	-0.8011485E-02	0.4575168E-03	-0.1095941E-04
			v= 14
			-0.2542487E+02	0.6270535E+01	-0.1908366E+01	0.1817239E+00	0.4856512E-01	-0.1807688E-01	0.2474441E-02	-0.1619281E-03	0.4199320E-05

.. toctree::

.. toctree::
    :hidden:
