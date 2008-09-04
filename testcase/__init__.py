__all__=['CAMxWind', 'CAMxLandUse', 'CAMxTemperature', 'CAMxVerticalDiffusivity', 'CAMxHumidity', 'CAMxHeightPressure', 'CAMxCloudRain', 'CAMxAverage', 'CAMxAreaEmissions', 'CAMxPointSource',]
from os.path import join

CAMxWind=join(__path__[0],'utils/CAMx/wind/camx_wind.20000825.hgbpa_04km.TCEQuh1_eta.v43')
CAMxLandUse=join(__path__[0],'utils/CAMx/landuse/camx_landuse.hgbpa_04km')
CAMxTemperature=join(__path__[0],'utils/CAMx/temperature/camx_temp.20000825.hgbpa_04km.TCEQuh1_eta.v43')
CAMxVerticalDiffusivity=join(__path__[0],'utils/CAMx/vertical_diffusivity/camx_kv.20000825.hgbpa_04km.TCEQuh1_eta.v43.tke')
CAMxHumidity=join(__path__[0],'utils/CAMx/humidity/camx_hum.20000825.hgbpa_04km.TCEQuh1_eta.v43')
CAMxHeightPressure=join(__path__[0],'utils/CAMx/height_pressure/camx_zp.20000825.hgbpa_04km.TCEQuh1_eta.v43')
CAMxCloudRain=join(__path__[0],'utils/CAMx/cloud_rain/camx_cr.20000825.hgbpa_04km.TCEQuh1_eta.v43')
CAMxAverage=join(__path__[0],'utils/CAMx/uamiv/camx420_cb4.20000816.hgb8h.base1b.psito2n2.TCEQuh1_eta_tke.avrg')
CAMxAreaEmissions=join(__path__[0],'utils/CAMx/uamiv/camx_cb4_ei_lo.20000825.hgb8h.base1b.psito2n2.hgbpa_04km')
CAMxPointSource=join(__path__[0],'utils/CAMx/point_source/camx_cb4_ei_el.20000825.hgb8h.base1b.psito2n2')