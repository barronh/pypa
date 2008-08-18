#!/usr/bin/python
""" colorscale.py version 0.11 2005-05-16 2005-05-22
Copyright Byeong-Uk Kim

Modifed by HEJ, 2006-09-14 to add a 'no_difference' value and color 
to colorscale for use with difference tile plots.
Default values are None and the operation is ignored.

This module provides a class for flexible color scale object.
Instance construction requires user-supplied tuple of tile intervals, 
  label positions, color tables.
As an option, user also can define color for values beyond the tile intervals.
For differece plot scales, user can define a "no difference" boundary
and a color to return when abs(differences) are less than this value.
"""

# python lib
from warnings import warn
from numarray import *
import numarray.ieeespecial as ieee # for inf, -inf
from string import *
from bisect import *

# external lib
from ChartDirector import *         # Chart Director


# begin colorscale

class colorscale:
    def __init__(self, color_code = (0x0000ff),  
                 steps = (0.0, 1.0), labels = (0.0, 1.0),
                 minus_inf_color = 0xffffff, plus_inf_color = 0x000000, 
                 label_fmt = "3.0f", 
                 diff_less_than_value = None, 
                 diff_less_than_color = None):
        # Setup whole color code schemes including two extreme value cases
        self.color_code = list(color_code)
        self.color_code.insert(0,minus_inf_color)
        self.color_code.append(plus_inf_color)
        
        self.min = steps[0]
        self.max = steps[-1]
        self.steps = list(steps)
        self.steps.insert(0,ieee.minus_inf)
        self.steps.append(ieee.plus_inf)
        
        self.diff_less_than_value = diff_less_than_value
        self.diff_less_than_color = diff_less_than_color
        self.labels = ['']*len(self.steps)
        self.labels[0] = ''
        self.labels[-1] = ''  
        if self.diff_less_than_value :   # make sure 0 val is 0.0
            self.steps[bisect(steps, diff_less_than_value)] = 0.0
        for label_idx in labels:
            self.labels[label_idx+1] = (" % "+label_fmt) % self.steps[label_idx+1]
    
    def tilecolor(self,val = 0.):
        if ( self.diff_less_than_value != None and self.diff_less_than_color != None):
            if abs(val) < self.diff_less_than_value :
                return self.diff_less_than_color
        for idx in xrange(0,len(self.steps)-1):
            if (self.steps[idx]<=val and val<self.steps[idx+1]):
                return self.color_code[idx]
    
    def colorbarplot(self, showinf = False):
        # return ChartDirector obj for adding color bars to various plots
        plot_width = 50.
        plot_height = 200.
        bar_width = 10.
        bar_height = 200.
        num_colors = len(self.color_code)
        each_color_height = ceil(bar_height/num_colors)
        self.colorbar = XYChart(plot_width, plot_height, 
                                0xffffff, 0xffffff, 0)
        self.colorbar.setYAxisOnRight(True) 
        self.colorbar.setPlotArea(0, 0, bar_width, bar_height, 
                                    Transparent, Transparent, Transparent, 
                                    Transparent, Transparent)
        self.colorbar.xAxis().setColors(Transparent, Transparent, 
                                    Transparent, Transparent)
        self.colorbar.xAxis().setLinearScale(0, bar_width)
        self.colorbar.yAxis().setColors(Transparent, Transparent, 
                                    Transparent, Transparent)
        self.colorbar.yAxis().setLabelStyle("normal", 8, 0x000000)
        self.colorbar.yAxis().setLinearScale(0, bar_height)
        self.colorbar.yAxis().setLabels(self.labels)
                
        if showinf == False:
            colorbaroffset = 1
        else:
            colorbaroffset = 0

        for idx in xrange(colorbaroffset, num_colors-colorbaroffset):
            each_color = self.color_code[idx]
            self.colorbar.addScatterLayer([bar_width],[idx+.5], "", SquareSymbol,
                    bar_width, each_color, each_color).setSymbolScale([bar_width],
                                    PixelScale,[each_color_height],PixelScale)
        
        return self.colorbar

# UNC's prefered color scheme for regular O3 tile plot
#  0.0 to 180.0 ppb in 
#    6 major steps (colors) white, green, yellow, orange, red, purple
#    of 5 shades of each color gives 30 ppb per color
#    with tics at 0, 30, 60, 90, 120, 150, 180
#  
color_table_30 = (0xffffff, 0xe6e6e6, 0xc8c8c8, 0xafafaf, 0x969696, 
                  0x00ff46, 0x00e646, 0x00d246, 0x00af46, 0x007d46, 
                  0xf5f500, 0xebeb00, 0xe1e100, 0xd7d700, 0xc8c800, 
                  0xffb400, 0xffa000, 0xff8c00, 0xff7800, 0xff6400, 
                  0xffbebe, 0xff9696, 0xff7878, 0xff5050, 0xff0000, 
                  0xff00ff, 0xdc00ff, 0xc800ff, 0xaf00ff, 0x9b00ff)

max_val = 300.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_30)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_30)+5, 5))
max300step6 = colorscale(color_table_30, steps, labels)

max_val = 180.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_30)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_30)+5, 5))
max180step6 = colorscale(color_table_30, steps, labels)

max_val = 120.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_30)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_30)+5, 5))
max120step6 = colorscale(color_table_30, steps, labels)

# Color scheme for O3 tile plots which can accomodate negative values
#  for difference plots.
color_table_28_negdiff = (
                  0x0000C0, 0x0000EE, 0x003EFF, 0x007FFF, # blues
                  0xe6e6e6, 0xcfcfcf, 0xafafaf, 0x969696, # greys
                  0x00ff46, 0x00e646, 0x00af46, 0x007d46, # greens
                  0xf5f500, 0xebeb00, 0xd7d700, 0xc8c800, # yellows
                  0xffb400, 0xffa000, 0xff7800, 0xff6400, # oranges
                  0xffbebe, 0xff9696, 0xff5050, 0xff0000, # reds
                  0xff00ff, 0xdc00ff, 0xaf00ff, 0x9b00ff  # purples
                  )

color_table_24_negdiff = (
                  0x0000C0, 0x0000EE, 0x003EFF, 0x007FFF, # blues
                  0xe6e6e6, 0xcfcfcf, 0xafafaf, 0x969696, # greys
                  0x00ff46, 0x00e646, 0x00af46, 0x007d46, # greens
                  0xf5f500, 0xebeb00, 0xd7d700, 0xc8c800, # yellows
                  0xffb400, 0xffa000, 0xff7800, 0xff6400, # oranges
                  0xffbebe, 0xff9696, 0xff5050, 0xff0000 # reds
                  )

color_table_30_negdiff = (
                  0x0000C0, 0x0000EE, 0x003EFF, 0x007FFF, 0x00CFFF, # blues
                  0xffffff, 0xe6e6e6, 0xc8c8c8, 0xafafaf, 0x969696, # greys
                  0x00ff46, 0x00e646, 0x00d246, 0x00af46, 0x007d46, # greens
                  0xf5f500, 0xebeb00, 0xe1e100, 0xd7d700, 0xc8c800, # yellows
                  0xffb400, 0xffa000, 0xff8c00, 0xff7800, 0xff6400, # oranges
                  0xffbebe, 0xff9696, 0xff7878, 0xff5050, 0xff0000) # purples
        
max_val = 300.
min_val = -50.
nodiff=2.
tilesize = (max_val-min_val)/len(color_table_28_negdiff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_28_negdiff)+4, 4))
diff_300maxstep6_50minstep1 = colorscale(color_table_28_negdiff, steps, labels,
                                 diff_less_than_value=nodiff,
                                 diff_less_than_color=0xebebeb, # colored 'gray' 
                                 minus_inf_color=0x000080)

max_val = 180.
min_val = -30.
nodiff=2.
tilesize = (max_val-min_val)/len(color_table_28_negdiff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_28_negdiff)+4, 4))
diff_180maxstep6_30minstep1 = colorscale(color_table_28_negdiff, steps, labels,
                                 diff_less_than_value=nodiff,
                                 diff_less_than_color=0xebebeb, # colored 'gray' 
                                 minus_inf_color=0x000080)

max_val = 120.
min_val = -20.
nodiff=1.
tilesize = (max_val-min_val)/len(color_table_28_negdiff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_28_negdiff)+4, 4))
diff_120maxstep6_20minstep1 = colorscale(color_table_28_negdiff, steps, labels,
                                 diff_less_than_value=nodiff,
                                 diff_less_than_color=0xebebeb, # colored 'gray' 
                                 minus_inf_color=0x000080)

max_val = 30.
min_val = -5.
nodiff  = 1.
tilesize = (max_val-min_val)/len(color_table_28_negdiff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_28_negdiff)+4, 4))
diff_30maxstep6_5minstep1_negblue = colorscale(color_table_28_negdiff, steps, labels,
                                               diff_less_than_value=nodiff,
                                               diff_less_than_color=0xebebeb, # colored 'gray' 
                                               minus_inf_color=0x000080)

max_val = 50.
min_val = -10.
nodiff  = 1.
tilesize = (max_val-min_val)/len(color_table_24_negdiff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_24_negdiff)+4, 4))
diff_50maxstep5_10minstep1_negblue = colorscale(color_table_24_negdiff, steps, labels,
                                               diff_less_than_value=nodiff,
                                               diff_less_than_color=0xebebeb, # colored 'gray' 
                                               minus_inf_color=0x000080)

max_val = 20.
min_val = -4.
nodiff  = 1.
tilesize = (max_val-min_val)/len(color_table_24_negdiff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_24_negdiff)+4, 4))
diff_20maxstep5_4minstep1_negblue = colorscale(color_table_24_negdiff, steps, labels,
                                               diff_less_than_value=nodiff,
                                               diff_less_than_color=0xebebeb, # colored 'gray' 
                                               minus_inf_color=0x000080)

max_val = 10.
min_val = -2.
nodiff  = 1.
tilesize = (max_val-min_val)/len(color_table_24_negdiff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_24_negdiff)+4, 4))
diff_10maxstep5_2minstep1_negblue = colorscale(color_table_24_negdiff, steps, labels,
                                               diff_less_than_value=nodiff,
                                               diff_less_than_color=0xebebeb, # colored 'gray' 
                                               minus_inf_color=0x000080)

color_table_25 = (0xffffff, 0xe6e6e6, 0xc8c8c8, 0xafafaf, 0x969696,
                  0x00ff46, 0x00e646, 0x00d246, 0x00af46, 0x007d46,
                  0xf5f500, 0xebeb00, 0xe1e100, 0xd7d700, 0xc8c800,
                  0xffb400, 0xffa000, 0xff8c00, 0xff7800, 0xff6400,
                  0xffbebe, 0xff9696, 0xff7878, 0xff5050, 0xff0000)
        
# UNC's prefered color scheme for max 100
max_val = 100.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max100step5 = colorscale(color_table_25, steps, labels)

# UNC's prefered color scheme for regular NOy and NOx tile plot
# 0.0 to 50.0 ppb in
#   5 major steps (colors) white, green, yellow, orange, red
#   of 5 shades of each color gives 10 ppb per color
#   with tics at 0, 10, 20, 30, 40, 50
#
max_val = 50.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max50step5 = colorscale(color_table_25, steps, labels)

#  UNC option for NOy NOx tile plots
max_val = 40.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max40step5 = colorscale(color_table_25, steps, labels)

nodiff = 1.0 # removes floating point errors +/- nodiff
max40step5_negblue = colorscale(color_table_25, steps, labels,
                                diff_less_than_value=nodiff,
                                diff_less_than_color=0xebebeb, # colored 'gray' 
                                minus_inf_color=0x000080)


#  UNC's prefered color scheme for regular FORM
max_val = 20.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max20step5 = colorscale(color_table_25, steps, labels)

#  UNC's prefered color scheme for regular ALD2
max_val = 15.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max15step5 = colorscale(color_table_25, steps, labels)

#  UNC's prefered color scheme for regular CO
max_val = 1000.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max1000step5 = colorscale(color_table_25, steps, labels)

#  UNC optional color scheme for regular CO
max_val = 800.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max800step5 = colorscale(color_table_25, steps, labels)

#  UNC's prefered color scheme for regular PAR
max_val = 500.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max500step5 = colorscale(color_table_25, steps, labels)

#  UNC optional color scheme for regular CO
max_val = 400.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max400step5 = colorscale(color_table_25, steps, labels)

#  UNC optional color scheme for ETH
max_val = 80.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max80step5 = colorscale(color_table_25, steps, labels)

#  UNC's prefered color scheme for regular ETH, OLE, and ISOP
max_val = 10.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max10step5 = colorscale(color_table_25, steps, labels)

#  UNC's prefered color scheme for trace gases
max_val = 1.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max1step5 = colorscale(color_table_25, steps, labels, label_fmt = "3.1f")

#  UNC's experimental color scheme for imputed HONO
max_val = 2.
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max2step5 = colorscale(color_table_25, steps, labels, label_fmt = "3.1f")

#  UNC's experimental color scheme for HONO
max_val = 0.1
min_val = 0.
tilesize = (max_val-min_val)/len(color_table_25)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25)+5, 5))
max01step5 = colorscale(color_table_25, steps, labels, label_fmt = "3.2f")

color_table_30_diff = (
        0x551a8b, 0x660198, 0x820bbb, 0xb23aee, 0xbf5fff,
        0x3232cd, 0x2e37fe, 0x4876ff, 0x5cacee, 0x74bbfb,
        0x007d46, 0x00af46, 0x00d246, 0x00e646, 0x00ff46, 
        0xf5f500, 0xebeb00, 0xe1e100, 0xd7d700, 0xc8c800, 
        0xffb400, 0xffa000, 0xff8c00, 0xff7800, 0xff6400,
        0xffbebe, 0xff9696, 0xff7878, 0xff5050, 0xff0000)

# UNC's prefered color scheme for 30 ppb differential tile plot
# -30 to 30 ppb in 
#   15 positive yellow, orange, red
#   15 negative purple, blue, green
#   with tics at -30, -20, -10, 0, 10, 20, 30
#   Set +1 to -1 ppb as a no-difference to be colored 'gray'

max_val = 30.
min_val = -30.
tilesize = (max_val-min_val)/len(color_table_30_diff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_30_diff)+5, 5))
nodiff = 1.0 # removes floating point errors +/- nodiff
diff30step6 = colorscale(color_table_30_diff, steps, labels, 
                         diff_less_than_value=nodiff,
                         diff_less_than_color=0xebebeb) # colored 'gray' 
diff30step6_negblue = colorscale(color_table_30_diff, steps, labels, 
                                 diff_less_than_value=nodiff,
                                 diff_less_than_color=0xebebeb, # colored 'gray' 
                                 minus_inf_color=0x000080)


# UNC's prefered color scheme for 10 ppb differential tile plot
# -10 to 10 ppb in 
#   15 positive yellow, orange, red
#   15 negative purple, blue, green
#   with tics at -10, -5, 0, 5, 10
#   Set +0.5 to -0.5 ppb as a no-difference to be colored 'gray'

max_val = 10.
min_val = -10.
tilesize = (max_val-min_val)/len(color_table_30_diff)
steps  = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_30_diff)+5, 5))
nodiff = 0.5 # removes floating point errors +/- nodiff
diff10step6 = colorscale(color_table_30_diff, steps, labels, 
                         label_fmt = "5.1f",
                         diff_less_than_value=nodiff,
                         diff_less_than_color=0xebebeb) # colored 'gray' 

diff10step6_negblue = colorscale(color_table_30_diff, steps, labels, 
                                 label_fmt = "5.1f",
                                 diff_less_than_value=nodiff,
                                 diff_less_than_color=0xebebeb, # colored 'gray' 
                                 minus_inf_color=0x000080)

color_table_30_diff = (
        0x551a8b, 0x660198, 0x820bbb, 0xb23aee, 0xbf5fff,
        0x3232cd, 0x2e37fe, 0x4876ff, 0x5cacee, 0x74bbfb,
        0x007d46, 0x00af46, 0x00d246, 0x00e646, 0x00ff46, 
        0xf5f500, 0xebeb00, 0xe1e100, 0xd7d700, 0xc8c800, 
        0xffb400, 0xffa000, 0xff8c00, 0xff7800, 0xff6400,
        0xffbebe, 0xff9696, 0xff7878, 0xff5050, 0xff0000)

max_val = 40.
min_val = -40.
tilesize = (max_val-min_val)/len(color_table_30_diff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_30_diff)+5, 5))
nodiff = 1.0 # removes floating point errors +/- nodiff
diff40step6 = colorscale(color_table_30_diff, steps, labels, 
                         diff_less_than_value=nodiff,
                         diff_less_than_color=0xebebeb) # colored 'gray' 
diff40step6_negblue = colorscale(color_table_30_diff, steps, labels, 
                                 diff_less_than_value=nodiff,
                                 diff_less_than_color=0xebebeb, # colored 'gray' 
                                 minus_inf_color=0x000080)

color_table_25_diff = (
                  0x551a8b, 0x660198, 0x820bbb, 0xb23aee, 0xbf5fff,
                  0x3232cd, 0x2e37fe, 0x4876ff, 0x5cacee, 0x74bbfb,
                  0x00d246, 0x00e646, 0xc8c8c8, 0xebeb00, 0xe1e100, 
                  0xffb400, 0xffa000, 0xff8c00, 0xff7800, 0xff6400,
                  0xffbebe, 0xff9696, 0xff7878, 0xff5050, 0xff0000)

max_val = 40.
min_val = -40.
tilesize = (max_val-min_val)/len(color_table_25_diff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25_diff)+5, 5))
nodiff = 1.0 # removes floating point errors +/- nodiff
diff40step5 = colorscale(color_table_25_diff, steps, labels, 
                         diff_less_than_value=nodiff,
                         diff_less_than_color=0xebebeb) # colored 'gray' 
diff40step5_negblue = colorscale(color_table_25_diff, steps, labels, 
                                 diff_less_than_value=nodiff,
                                 diff_less_than_color=0xebebeb, # colored 'gray' 
                                 minus_inf_color=0x000080)

max_val = 20.
min_val = -20.
tilesize = (max_val-min_val)/len(color_table_25_diff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25_diff)+5, 5))
nodiff = 0.25
diff20step5 = colorscale(color_table_25_diff, steps, labels,
                         diff_less_than_value=nodiff,
                         diff_less_than_color=0xebebeb) # colored 'gray' 
diff20step5_negblue = colorscale(color_table_25_diff, steps, labels, 
                                 diff_less_than_value=nodiff,
                                 diff_less_than_color=0xebebeb, # colored 'gray' 
                                 minus_inf_color=0x000080)

max_val = 10.
min_val = -10.
tilesize = (max_val-min_val)/len(color_table_25_diff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25_diff)+5, 5))
nodiff = 0.1
diff10step5 = colorscale(color_table_25_diff, steps, labels,
                         diff_less_than_value=nodiff,
                         diff_less_than_color=0xebebeb) # colored 'gray' 
diff10step5_negblue = colorscale(color_table_25_diff, steps, labels, 
                                 diff_less_than_value=nodiff,
                                 diff_less_than_color=0xebebeb, # colored 'gray' 
                                 minus_inf_color=0x000080)
                                 
                                 
max_val = 0.5
min_val = -0.5
tilesize = (max_val-min_val)/len(color_table_25_diff)
steps = tuple(arange(min_val, max_val+tilesize, tilesize))
labels = tuple(arange(0, len(color_table_25_diff)+5, 5))
nodiff = 0.01
diff05step5 = colorscale(color_table_25_diff, steps, labels, label_fmt = "3.2f",
                         diff_less_than_value=nodiff,
                         diff_less_than_color=0xebebeb) # colored 'gray' 
diff05step5_negblue = colorscale(color_table_25_diff, steps, labels, label_fmt = "3.2f",
                                 diff_less_than_value=nodiff,
                                 diff_less_than_color=0xebebeb, # colored 'gray' 
                                 minus_inf_color=0x000080)
# Cleaning up local vars
major_tilesize = None
minor_tilesize = None
max_val = None
min_val = None
ntile = None
nlabel = None
ncolor = None

# Ditionary for the predefined color scales
ColorScaleDict={ "max1000step5": max1000step5 ,
                 "max800step5" : max800step5  ,
                 "max500step5" : max500step5  ,
                 "max400step5" : max400step5  ,
                 "max300step6" : max300step6  ,
                 "max180step6" : max180step6  ,
                 "max120step6" : max120step6  ,
                 "max80step5"  : max80step5  ,
                 "max50step5"  : max50step5   ,
                 "max40step5"  : max40step5   ,
         "max40step5_negblue"  : max40step5_negblue,
                 "max20step5"  : max20step5   ,
                 "max15step5"  : max15step5   ,
                 "max10step5"  : max10step5   ,
                 "max2step5"   : max2step5    ,
                 "max1step5"   : max1step5    ,
                 "max01step5"  : max01step5   ,
                 "diff30step6" : diff30step6  ,
                 "diff10step6" : diff10step6  ,
         "diff10step6_negblue" : diff10step6_negblue,
                 "diff40step5" : diff40step5  ,
                 "diff20step5" : diff20step5  ,
                 "diff10step5" : diff10step5  ,
                 "diff05step5" : diff05step5  ,
         "diff05step5_negblue" : diff05step5_negblue,
         "diff10step5_negblue" : diff10step5_negblue,
         "diff20step5_negblue" : diff20step5_negblue,
         "diff40step5_negblue" : diff40step5_negblue,
 "diff_120maxstep6_20minstep1" : diff_120maxstep6_20minstep1,
 "diff_180maxstep6_30minstep1" : diff_180maxstep6_30minstep1,
 "diff_300maxstep6_50minstep1" : diff_300maxstep6_50minstep1,
"diff_30maxstep6_5minstep1_negblue" : diff_30maxstep6_5minstep1_negblue,
"diff_50maxstep5_10minstep1_negblue" : diff_50maxstep5_10minstep1_negblue,
"diff_20maxstep5_4minstep1_negblue" : diff_20maxstep5_4minstep1_negblue,
"diff_10maxstep5_2minstep1_negblue" : diff_10maxstep5_2minstep1_negblue}
warn('ColorScale is a "borrowed" file.  When pyPASS redevelopment is complete, this should be imported rather than copied')
# Test code runs when the module called by its name
if __name__ == '__main__':
    print "Test code is running..."
    test_conc = 5.
    print " Value of %4.1f corresponds to RGB code of %#x by max180step6" % (test_conc, max180step6.tilecolor(test_conc))
    print " Value of %4.1f corresponds to RGB code of %#x by max50step5" % (test_conc, max50step5.tilecolor(test_conc))
    print " Value of %4.1f corresponds to RGB code of %#x by diff30step6" % (test_conc, diff30step6.tilecolor(test_conc))

    for label, scale in ColorScaleDict.items() :
        scale.colorbarplot().makeChart(label+".png")
        if "_neg" in label : scale.colorbarplot(showinf=True).makeChart(label+".png")
