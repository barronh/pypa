#!/usr/bin/env python
__all__ = ['phy_plot']

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from net_balance import get_pure_mech
from datetime import datetime
from pyPA.utils.util import AttrDict
from pyPA.netcdf import NetCDFFile
from numpy import array, concatenate, zeros
from pylab import figure, title, plot_date, savefig, legend, axis, grid, axes, xlabel, ylabel
from matplotlib.dates import DateFormatter
from matplotlib.colors import cnames
from matplotlib.font_manager import FontProperties
import re

def get_date(mrg_file, job):
    date_ints = mrg_file.variables['TFLAG'][:,0,:]
    date_objs = array([datetime.strptime("%iT%06i" % (d,t), "%Y%jT%H%M%S") for d,t in date_ints])
    if job.end_date:
        date_objs = concatenate([date_objs[[0]]-(date_objs[-1]-date_objs[-2]), date_objs])
    else:
        date_objs = concatenate([date_objs, date_objs[[-1]]+(date_objs[-1]-date_objs[-2])])

    date_objs = date_objs.repeat(2,0)[1:-1]
    return date_objs

def chem_plot(job, nlines = 8, fmt = 'pdf'):
    mech = get_pure_mech('_'.join([job.mechanism,job.model]))
    mrg_file = NetCDFFile(job.mrgfile,'r')
    mech.set_mrg(mrg_file)
    
    units = mrg_file.variables['IRR'].units
    date_objs = get_date(mrg_file, job)

    for species in job.species:
        fig = figure()
        grid(True)
        title(job.title % locals())
                
        for species in job.species:
            chem_colors = ['indigo', 'gold', 'firebrick', 'indianred', 'yellow', 'darkolivegreen', 'darkseagreen', 'slategrey', 'darkslategrey', 'mediumvioletred', 'mediumorchid', 'chartreuse', 'mediumslateblue', 'black', 'springgreen', 'crimson', 'lightsalmon', 'brown', 'turquoise', 'olivedrab', 'cyan', 'silver', 'skyblue', 'gray', 'darkturquoise', 'goldenrod', 'darkgreen', 'darkviolet', 'darkgray', 'lightpink', 'teal', 'darkmagenta', 'lightgoldenrodyellow', 'lavender', 'yellowgreen', 'thistle', 'violet', 'navy', 'dimgrey', 'orchid', 'blue', 'ghostwhite', 'honeydew', 'cornflowerblue', 'darkblue', 'darkkhaki', 'mediumpurple', 'cornsilk', 'red', 'bisque', 'slategray', 'darkcyan', 'khaki', 'wheat', 'deepskyblue', 'darkred', 'steelblue', 'aliceblue', 'lightslategrey', 'gainsboro', 'mediumturquoise', 'floralwhite', 'coral', 'purple', 'lightgrey', 'lightcyan', 'darksalmon', 'beige', 'azure', 'lightsteelblue', 'oldlace', 'greenyellow', 'royalblue', 'lightseagreen', 'mistyrose', 'sienna', 'lightcoral', 'orangered', 'navajowhite', 'lime', 'palegreen', 'burlywood', 'seashell', 'mediumspringgreen', 'fuchsia', 'papayawhip', 'blanchedalmond', 'peru', 'aquamarine', 'white', 'darkslategray', 'ivory', 'dodgerblue', 'lemonchiffon', 'chocolate', 'orange', 'forestgreen', 'slateblue', 'olive', 'mintcream', 'antiquewhite', 'darkorange', 'cadetblue', 'moccasin', 'limegreen', 'saddlebrown', 'grey', 'darkslateblue', 'lightskyblue', 'deeppink', 'plum', 'aqua', 'darkgoldenrod', 'maroon', 'sandybrown', 'magenta', 'tan', 'rosybrown', 'pink', 'lightblue', 'palevioletred', 'mediumseagreen', 'dimgray', 'powderblue', 'seagreen', 'snow', 'mediumblue', 'midnightblue', 'palegoldenrod', 'whitesmoke', 'darkorchid', 'salmon', 'lightslategray', 'lawngreen', 'lightgreen', 'tomato', 'hotpink', 'mediumaquamarine', 'green', 'blueviolet', 'darkgrey']

            reactions = mech.find_rxns(mech(species), mech(species), False)
            reactions = [ (abs(mech('%s[%s]' % (rxn, species))).sum(),rxn) for rxn in reactions]
            reactions.sort(reverse = True)
            reactions = [r for v,r in reactions]
            
            options = {}
            options.setdefault('ls', '-')
            options.setdefault('lw', 3)
            options.setdefault('marker', 'None')

            for rxn in reactions[:nlines-1]:
                data = mech('%s[%s]' % (rxn, species))
                options['color'] = chem_colors.pop()
                options['label'] = re.compile('\d.\d{5}\*').sub('', str(mech.reaction_dict[rxn]))
                plot_date(date_objs, data.repeat(2,0), **options)

            data = zeros(data.shape, dtype = data.dtype)
            for rxn in reactions[nlines-1:]:
                data += mech('%s[%s]' % (rxn, species))
            options['color'] = chem_colors.pop()
            options['label'] = 'other'
            plot_date(date_objs, data.repeat(2,0), **options)
            
            options['color'] = 'black'
            options['marker'] = 'o'
            options['label'] = 'Chem'
            data = mech('%s[%s]'  % (job.chem, species)).array().repeat(2,0)
            plot_date(date_objs, data, **options)

            ax = axes()
            xlabel('Time')
            ylabel(units)
            ax.xaxis.set_major_formatter(DateFormatter('%H'))
            fig.autofmt_xdate()
            legend(loc='lower right', prop = FontProperties(size=10))
            savefig('%s_IRR.%s' % (species, fmt), format = fmt)
    
def phy_plot(job, filter = True, fmt = 'pdf'):
    mech = get_pure_mech('_'.join([job.mechanism,job.model]))
    mrg_file = NetCDFFile(job.mrgfile,'r')
    mech.set_mrg(mrg_file)
    
    units = mrg_file.variables['IPR'].units
    date_objs = get_date(mrg_file, job)
    for species in job.species:
        fig = figure()
        grid(True)
        title(job.title % locals())
        options = {}
        options.setdefault('c', 'k')
        options.setdefault('ls', '-')
        options.setdefault('lw', 3)
        options.setdefault('marker', 'o')
        data = mech('%s[%s]'  % (job.init, species)).array()
        plot_date(date_objs[::2], data, **options)
        options.setdefault('label', 'Conc')
        options['marker'] = 'x'
        data = mech('%s[%s]'  % (job.final, species)).array()
        plot_date(date_objs[1::2], data, **options)
        
        for process in job.process.keys():
            options = job.process[process]
            options.setdefault('linestyle', '-')
            options.setdefault('linewidth', 3)
            options.setdefault('label', process)
            options.setdefault('marker', 'None')
            options.setdefault('color', 'None')
            data = mech('%s[%s]' % (process, species)).array().repeat(2,0)
            if data.nonzero()[0].any() or not filter:
                plot_date(date_objs, data, **options)
        options = job.species[species]
        axis(**options)
        ax = axes()
        xlabel('Time')
        ylabel(units)
        ax.xaxis.set_major_formatter(DateFormatter('%H'))
        fig.autofmt_xdate()
        legend()
        savefig('%s_IPR.%s' % (species, fmt), format = fmt)
        
if __name__ == '__main__':
    from pyPA.graphing.phy_plot import phy_plot
    job = AttrDict()
    job.mechanism = 'cb05'
    job.model = 'camx'
    job.mrgfile = '/Users/barronh/Development/net_balance/src/net_balance/testdata/test.mrg.nc'
    job.species  = dict(NO2=dict(ymin=None))
    job.process = dict(H_Trans=dict(color = 'b'), 
                       V_Trans=dict(color = 'g'),
                       Emissions=dict(color = 'r'),
                       Deposit=dict(color = 'c'),
                       Aero_Chem=dict(color = 'm'),
                       CHEM=dict(color = 'y'),
                       Motion=dict(color = 'teal'))
    job.init = 'INIT'
    job.final = 'FCONC'
    job.chem = 'CHEM'
    job.title = '%(species)s IPR plot'
    job.end_date = False
    phy_plot(job)
    chem_plot(job)