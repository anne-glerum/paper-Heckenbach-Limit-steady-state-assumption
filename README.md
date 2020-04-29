# paper-Heckenbach-Limit-steady-state-assumption

This repository belongs to the paper

    Is there a speed limit for the thermal steady-state assumption
    in continental rifts?
    by
    E. L Heckenbach, S. Brune, A. C. Glerum and J. Bott.

and contains an input file and code to reproduce the computations in the paper.

Contents
--------
``aspect_plugins_time_dependent/``

The ASPECT plugins created specifically for the time-dependent models of this paper. They should be build as shared libraries in conjunction with ASPECT 2.0.0-pre commit 585d1c3c99de259057408ea90aab5dbe963ecb40, see the supplied installation instructions. A full ASPECT branch including these plugins is available at ``https://github.com/anne-glerum/aspect/tree/paper-Heckenbach-Limit-steady-state-assumption``.

``aspect_plugins_steady_state/``

The ASPECT plugins created specifically for the steady-state models of this paper. They should be build as shared libraries in conjunction with ASPECT 2.0.0-pre commit 791f903229e4cdc65c04710fd0d9211a7250948d, see the supplied installation instructions. A full ASPECT branch including these plugins is available at ``https://github.com/anne-glerum/aspect/tree/paper-Heckenbach-Limit-steady-state-assumption-2``.

``input_file_time_dependent.prm``

The commented ASPECT input file used to create a representative computation in the paper (extension velocity of 4 mm/yr and a 30 km thick crust). 

``input_file_steady_state.prm``

The commented ASPECT input file used to create a representative steady-state computation in the paper. The required input data files were created from the time dependent model run with ``input_file_time_dependent.prm``. 

``input_steady_state/``
The folder containing the input data files for the steady-state computation, including surface topography and upper crust/lower crust, Moho and LAB interface depths. 

``ASPECT_install_configuration.txt``

The specific installation configurations for ASPECT and deal.ii as used for the paper.
