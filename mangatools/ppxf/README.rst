
**pPXF: Full Spectrum Fitting of Galactic and Stellar Spectra**

.. image:: https://img.shields.io/pypi/v/ppxf.svg
    :target: https://pypi.org/project/ppxf/
.. image:: https://img.shields.io/badge/arXiv-1607.08538-orange.svg
    :target: https://arxiv.org/abs/1607.08538
.. image:: https://img.shields.io/badge/DOI-10.1093/mnras/stw3020-green.svg
        :target: https://doi.org/10.1093/mnras/stw3020

This ``pPXF`` package contains a Python implementation of the Penalized
PiXel-Fitting (``pPXF``) method to fit the stellar and gas kinematics,
as well as the stellar population of galaxies. The method was originally
described in `Cappellari & Emsellem (2004)
<http://adsabs.harvard.edu/abs/2004PASP..116..138C>`_
and was significantly upgraded in `Cappellari (2017)
<http://adsabs.harvard.edu/abs/2017MNRAS.466..798C>`_.

.. contents:: :depth: 1

Attribution
-----------

If you use this software for your research, please cite at least
`Cappellari (2017) <http://adsabs.harvard.edu/abs/2017MNRAS.466..798C>`_.
The BibTeX entry for the paper is::

    @ARTICLE{Cappellari2017,
        author = {{Cappellari}, M.},
        title = "{Improving the full spectrum fitting method:
            accurate convolution with Gauss-Hermite functions}",
        journal = {MNRAS},
        eprint = {1607.08538},
        year = 2017,
        volume = 466,
        pages = {798-811},
        doi = {10.1093/mnras/stw3020}
    }

Installation
------------

install with::

    pip install ppxf

Without write access to the global ``site-packages`` directory, use::

    pip install --user ppxf

Usage Examples
--------------

To learn how to use the main program ``pPXF`` run the example programs in the
``ppxf/examples`` directory and read the detailed documentation in the file
``ppxf.py`` or as PDF from `<http://purl.org/cappellari/software>`_.

