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

def run():
	unittest.main()

if __name__=='__main__':
	run()

