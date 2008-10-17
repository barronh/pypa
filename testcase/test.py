"""
Need to add docs
"""
import unittest
from pyPA.utils.CAMx.wind.Memmap import TestMemmap as WindMemmap
from pyPA.utils.CAMx.humidity.Memmap import TestMemmap as HumidityMemmap
from pyPA.utils.CAMx.temperature.Memmap import TestMemmap as TemperatureMemmap
from pyPA.utils.CAMx.vertical_diffusivity.Memmap import TestMemmap as VerticalDiffusivityMemmap
from pyPA.utils.CAMx.height_pressure.Memmap import TestMemmap as HeightPressureMemmap
from pyPA.utils.CAMx.landuse.Memmap import TestMemmap as LandUseMemmap
from pyPA.utils.CAMx.cloud_rain.Memmap import TestMemmap as CloudRainMemmap
from pyPA.utils.CAMx.uamiv.Memmap import TestMemmap as UAMIVMemmap
from pyPA.utils.CAMx.point_source.Memmap import TestMemmap as PointSourceMemmap
from pyPA.utils.CAMx.ipr.Memmap import TestMemmap as IPRMemmap
from pyPA.utils.CAMx.irr.Memmap import TestMemmap as IRRMemmap
from pyPA.pappt.kvextract import TestPbl
from pyPA.pappt.lagrangian import TestBoxes
from pyPA.pappt.pappt import TestExtractor
from pyPA.utils.ArrayTransforms import TestInteriorVertex
from pyPA.utils.FortranFileUtil import TestFileUtils
from pyPA.utils.sci_var import PseudoNetCDFTest
from pyPA.utils.util import CompareTime

def run():
	import pyPA.testcase.test
	suite=unittest.TestLoader().loadTestsFromModule(pyPA.testcase.test)
	unittest.TextTestRunner(verbosity=2).run(suite)

if __name__=='__main__':
	run()

