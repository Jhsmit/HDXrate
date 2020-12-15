=======
HDXrate
=======


.. image:: https://img.shields.io/pypi/v/hdxrate.svg
        :target: https://pypi.python.org/pypi/hdxrate

.. image:: https://img.shields.io/travis/Jhsmit/hdxrate.svg
        :target: https://travis-ci.com/Jhsmit/hdxrate

.. image:: https://readthedocs.org/projects/hdxrate/badge/?version=latest
        :target: https://hdxrate.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Python package collection for HDX intrinsic exchange rate calculation. 
The calculations are based on the following papers:

  Bai, Y., Milne, J. S., Mayne, L. & Englander, S. W. Primary structure effects on peptide group hydrogen exchange. `Proteins Structure, Function, and Bioinformatics <https://doi.org/10.1002/prot.340170110>`__ 17, 75–86 (1993)

  Connelly, G. P., Bai, Y., Jeng, M.-F. & Englander, S. W. Isotope effects in peptide group hydrogen exchange. `Proteins <http://doi.wiley.com/10.1002/prot.340170111>`__ 17, 87–92 (1993).


  Mori, S., Zijl, P. C. M. van & Shortle, D. Measurement of water–amide proton exchange rates in the denatured state of staphylococcal nuclease by a magnetization transfer technique. `Proteins Structure, Function, and Bioinformatics <https://doi.org/10.1002/(SICI)1097-0134(199707)28:3%3C325::AID-PROT3%3E3.0.CO;2-B>`__ 28, 325–332 (1997)
  
  Nguyen, D., Mayne, L., Phillips, M. C. & Walter Englander, S. Reference Parameters for Protein Hydrogen Exchange Rates. `J. Am. Soc. Mass Spectrom.  <https://pubs.acs.org/doi/abs/10.1021/jasms.8b05911>`__ 29, 1936–1939 (2018).

..
   1.Rubinson, K. A. Practical corrections for p(H,D) measurements in mixed H 2 O/D 2 O biological buffers. Anal. Methods 9, 2744–2750 (2017).


See also the excel sheet on the Englander group website: http://hx2.med.upenn.edu/download.html


* Free software: GNU General Public License v3


Features
--------

Calculate intrinsic rate of amide hydrogen exchange in proteins.

Installation
------------

::

   $ pip install hdxrate

or

::

   $ conda install -c conda-forge hdxrate

Usage
-----

::

   >>> from hdxrate import k_int_from_sequence
   >>> k_int_from_sequence('HHHHH', 300, 7.)
   array([0.00000000e+00, 2.62430718e+03, 6.29527446e+01, 6.29527446e+01,
       9.97734191e-01])


Credits
-------

HDXrate implementation is based on PSX intrinsic rate calculation
https://github.com/Niels-Bohr-Institute-XNS-StructBiophys/PSX

 Pedersen, M. C. et al. PSX, Protein–Solvent Exchange: software for calculation of deuterium-exchange effects in small-angle neutron scattering measurements from protein coordinates. `J Appl Cryst <https://doi.org/10.1107/S1600576719012469/>`__ 52, 1427–1436 (2019).



Maintenance
```````````

* Jochem Smit <jhsmit@gmail.com> / <jochem.smit@kuleuven.be>
