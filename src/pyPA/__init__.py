__doc__ = r"""
.. _pyPA
:mod:`pyPA` -- Python-based Process Analysis
============================================

.. module:: pyPA
   :platform: Unix, Windows
   :synopsis: Provides tools for analyzing Air Quality Model Process Analysis 
   data
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['utils','pappt', 'test']
if __name__ != '__main__':
    import utils
    import pappt
    import cmaq
    
    from test import run as test
else:
    from main import run
    run()
