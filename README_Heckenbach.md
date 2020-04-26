# paper-Heckenbach-Limit-steady-state-assumption

This repository belongs to the paper

    Is there a speed limit for the thermal steady-state assumption
    in continental rifts?
    by
    E. L Heckenbach, S. Brune, A. C. Glerum and J. Bott.

and contains an input file and code to reproduce the computations in the paper.

Contents
--------
``aspect_plugins_dynamic``

The ASPECT plugins created specifically for the time-dependent models of this paper. They should be build as shared libraries in conjunction with ASPECT 2.0.0-pre commit 585d1c3c99de259057408ea90aab5dbe963ecb40, see the supplied installation instructions. A full ASPECT branch including these plugins is available at ``https://github.com/anne-glerum/aspect/tree/paper-Heckenbach-Limit-steady-state-assumption``.

``aspect_plugins_steady_state``

The ASPECT plugins created specifically for the steady-state models of this paper. They should be build as shared libraries in conjunction with ASPECT 2.0.0-pre commit 791f903229e4cdc65c04710fd0d9211a7250948d, see the supplied installation instructions. A full ASPECT branch including these plugins is available at ``https://github.com/anne-glerum/aspect/tree/paper-Heckenbach-Limit-steady-state-assumption-2``.

``input_file.prm``

The commented ASPECT input file used to create the main computation in the paper. Please set the parameter ``Additional shared libaries`` to the library created from the supplied ASPECT plugins. 

``ASPECT_install_configuration.txt``

The specific installation configurations for ASPECT and deal.ii as used for the paper.
