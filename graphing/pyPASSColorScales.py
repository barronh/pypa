"""
This module is intended to make pyPASS color scales 
available for pyPA plotting.  pyPA plotting is mostly 
based on matplotlib, but pyPASS plotting is based on 
ChartDirector.  The colorscales have slightly different
formats.
"""
__all__=['DiscreteColormap','GradientColormap']
from ColorScale import ColorScaleDict as pyPASSColorScaleDict
from matplotlib.colors import ListedColormap, hex2color, LinearSegmentedColormap

def DiscreteColormap(name):
    out=ListedColormap(['#%06x' % i for i in pyPASSColorScaleDict[name].color_code][1:-1])
    out.set_under('#%06x' % pyPASSColorScaleDict[name].color_code[0])
    out.set_over('#%06x' % pyPASSColorScaleDict[name].color_code[-1])
    return out

def GradientColormap(name):
    colors=[hex2color(i) for i in ['#%06x' % i for i in pyPASSColorScaleDict[name].color_code][1:-1]]
    
    colord=dict(red=(),green=(),blue=())
    for ci,c in enumerate(['red','green','blue']):
        for index,bottom in enumerate(range(0,len(colors))):
            bottom/=float(len(colors))
            colord[c]=colord[c]+((bottom,colors[index][ci],colors[index][ci]),)
        
        colord[c]=colord[c]+((1,colors[index][ci],colors[index][ci]),)
    
    out=LinearSegmentedColormap(name,colord)
    out.set_under('#%06x' % pyPASSColorScaleDict[name].color_code[0])
    out.set_over('#%06x' % pyPASSColorScaleDict[name].color_code[-1])
    return out

for k in pyPASSColorScaleDict.keys():
    locals()[k+"_seg"]=DiscreteColormap(k)
    locals()[k]=GradientColormap(k)
    __all__.append(k)
    __all__.append(k+"_seg")