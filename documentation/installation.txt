Installation
============

This version of the code works on ubuntu 9.10. To install this code, download
pyphasesym-1.0, untar it and execute following commands on your terminal::

$ sudo python setup.py install

This allows you to use pyphasesym package as ::

import pyphasesym
or
from pyphasesym import *


Unittesting
-------------

In this package, we provide unittesting and regression testing. For 
unittesting exectute following command::

$ nosetests test_pyphasesym_unit.py

Regression Testing
--------------------

To carry out regression tests, you need matfiles with set of parameters
used on original matlab code of Dr. Peter Kovesi. Obtain his code from his
website. test_phasesym.m file is provided to generate testing files in matlab
to carry out regression testing. To generate test matfiles, open a matlab 
terminal (change to appropriate directory) and then execute following command 
on matlab prompt:

$ test_phasesym(n)

where n is number of iterations tests you want to carry out. The test_phasesym.m 
file selects the parameters at random for a given range of parameters for 
original phasesym code and writes matfiles. Store all the matlfiles in folder
matfiles/ directory of pyphasesym. Some example test matfiles are provided with
the package

Some matfiles for various parameters are stored in matfiles/ directory. Nosetests
can be used for regression tests as follows::

$ nosetests test_pyphasesym_regression.py


