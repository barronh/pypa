r"""
.. _pyPA
:mod:`pyPA` -- Python-based Process Analysis
============================================

.. module:: pyPA
   :platform: Unix, Windows
   :synopsis: Provides tools for analyzing Air Quality Model Process Analysis 
   data
.. moduleauthor:: Barron Henderson <barronh@unc.edu>
"""
__all__=['utils','pappt','graphing']

#from pyPA.pappt.loader import LoadPyPAFromYAML as run, LoadPAQAFromYAML as qa, template
try:
    from pyPA.testcase.test import run as test
except:
    def test():
        raise ValueError, "Testcase was not installed"

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_usage("Usage: %prog [-tq] [-m <model>] [-c <mechanism>]\n"+(" "*16)+" [-i <init name>] [-f <final name>] <yamlfile>")
    parser.add_option("-q", "--qa", dest="qa",action="store_true",default=False,
                        help="Check PA data", metavar="QAPA")
    
    parser.add_option("-t", "--template", dest="template",action="store_true",default=False,
                        help="Output template on standard out (configurable with -m and -c", metavar="Template")
    
    parser.add_option("-m", "--model", dest="model",default="new",
                        help="Model can either be camx, cmaq or wrfchem (for use with -t)", metavar="MODEL")
    
    parser.add_option("-c", "--mechanism", dest="mechanism", default="template",
                        help="Chemical mechanisms: cbivmech3 (for use with -t)", metavar="MECHANISM")
    
    parser.add_option("-i", "--initial", dest="initial", default="INIT")
    parser.add_option("-f", "--final", dest="final", default="FCONC")
    
    (options, args) = parser.parse_args()
    if options.template:
        from pyPA.pappt.loader import template
        template(options.model,options.mechanism)
        parser.exit()
    if len(args)<1:
        parser.error(msg="Requires a yaml file as an argument.  For a template use the -t option.  The template will be output to the stdout.")
    else:
        yamlpath=args[0]
    
    if len(args)==2:
        netbalance=args[1]
    
    if options.qa:
        qa(yamlpath)
    else:
        run(yamlpath)
