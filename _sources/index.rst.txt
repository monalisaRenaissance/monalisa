.. doc_renaissance documentation master file, created by
   sphinx-quickstart on Fri Sep 26 13:49:57 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================
MRI Reconstructions 
===================

From the Theory to the Implementation
=====================================

.. image:: ./images/Bastien.jpeg
   :width: 100 %


The `Monalisa <https://github.com/monalisaRenaissance/monalisa.git>`_ toolbox for MRI reconstruction has been originaly developed at CIBM-CHUV between 2018 and 2023 
by Bastien Milani.  

Originally, the development of the toolbox began by the implementation of non-cartesian reconstructions. 
The first reconstruction implemented was a gridded reconstruction which is part of the static non-iterative 
familly. After that, some static iterative reconstruction were added and later 3D-CINE iterative reconstructions
with temporal regularisation were implemented (4D and 5D), all for non-cartesian data. 
Iterative 3D-CINE reconstruction for cartesian data were then implemented on the same model.  

This toolbox can be used freely for any resonable personal or academic application, including publications. But any redistribution 
or commercial usage must be done in collaboration with lincensors (see the LICENSE file).  



.. toctree::
   :maxdepth: 3

   1_instalation
   2_introduction
   3_demos
   4_thermodynInfo





.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
