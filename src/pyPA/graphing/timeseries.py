#!/usr/bin/env python
__all__ = ['phy_plot']

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

from warnings import warn
from net_balance import get_pure_mech
from datetime import datetime
from pyPA.utils.util import AttrDict
from pyPA.netcdf import NetCDFFile
from numpy import array, concatenate, zeros
from pylab import figure, title, plot_date, savefig, legend, axis, grid, axes, xlabel, ylabel, subplot
from matplotlib.dates import DateFormatter
from matplotlib.colors import cnames
from matplotlib.font_manager import FontProperties
import re
import operator
import os

def get_date(mrg_file, conf):
    date_ints = mrg_file.variables['TFLAG'][:,0,:]
    date_objs = array([datetime.strptime("%iT%06i" % (d,t), "%Y%jT%H%M%S") for d,t in date_ints])
    if conf.end_date:
        date_objs = concatenate([date_objs[[0]]-(date_objs[-1]-date_objs[-2]), date_objs])
    else:
        date_objs = concatenate([date_objs, date_objs[[-1]]+(date_objs[-1]-date_objs[-2])])

    date_objs = date_objs.repeat(2,0)[1:-1]
    return date_objs

def chem_plot(conf, mech, date_objs, species, species_options, nlines, combine):
    units = mech.irr.units

    fig = figure()
    ax = axes([0.1,0.4,.8,.5])
    grid(True)
    title(conf.title % locals())
    chem_colors = ['indigo', 'gold', 'firebrick', 'indianred', 'yellow', 'darkolivegreen', 'darkseagreen', 'slategrey', 'darkslategrey', 'mediumvioletred', 'mediumorchid', 'chartreuse', 'mediumslateblue', 'black', 'springgreen', 'crimson', 'lightsalmon', 'brown', 'turquoise', 'olivedrab', 'cyan', 'silver', 'skyblue', 'gray', 'darkturquoise', 'goldenrod', 'darkgreen', 'darkviolet', 'darkgray', 'lightpink', 'teal', 'darkmagenta', 'lightgoldenrodyellow', 'lavender', 'yellowgreen', 'thistle', 'violet', 'navy', 'dimgrey', 'orchid', 'blue', 'ghostwhite', 'honeydew', 'cornflowerblue', 'darkblue', 'darkkhaki', 'mediumpurple', 'cornsilk', 'red', 'bisque', 'slategray', 'darkcyan', 'khaki', 'wheat', 'deepskyblue', 'darkred', 'steelblue', 'aliceblue', 'lightslategrey', 'gainsboro', 'mediumturquoise', 'floralwhite', 'coral', 'purple', 'lightgrey', 'lightcyan', 'darksalmon', 'beige', 'azure', 'lightsteelblue', 'oldlace', 'greenyellow', 'royalblue', 'lightseagreen', 'mistyrose', 'sienna', 'lightcoral', 'orangered', 'navajowhite', 'lime', 'palegreen', 'burlywood', 'seashell', 'mediumspringgreen', 'fuchsia', 'papayawhip', 'blanchedalmond', 'peru', 'aquamarine', 'white', 'darkslategray', 'ivory', 'dodgerblue', 'lemonchiffon', 'chocolate', 'orange', 'forestgreen', 'slateblue', 'olive', 'mintcream', 'antiquewhite', 'darkorange', 'cadetblue', 'moccasin', 'limegreen', 'saddlebrown', 'grey', 'darkslateblue', 'lightskyblue', 'deeppink', 'plum', 'aqua', 'darkgoldenrod', 'maroon', 'sandybrown', 'magenta', 'tan', 'rosybrown', 'pink', 'lightblue', 'palevioletred', 'mediumseagreen', 'dimgray', 'powderblue', 'seagreen', 'snow', 'mediumblue', 'midnightblue', 'palegoldenrod', 'whitesmoke', 'darkorchid', 'salmon', 'lightslategray', 'lawngreen', 'lightgreen', 'tomato', 'hotpink', 'mediumaquamarine', 'green', 'blueviolet', 'darkgrey']
    species_obj = mech(species)
    reactions = mech.find_rxns(species_obj, species_obj, False)
    reactions = [ rxn for rxn in reactions if rxn not in reduce(operator.add, combine) ]
    if combine != [()]:
        reactions = reactions + map(lambda t2: '+'.join(t2), combine)
    reactions = [ (abs(mech('(%s)[%s]' % (rxn, species))).sum(),rxn) for rxn in reactions]
    reactions.sort(reverse = True)
    reactions = [r for v,r in reactions]
    
    options = {}
    options.setdefault('ls', '-')
    options.setdefault('lw', 3)
    options.setdefault('marker', 'None')

    for rxn in reactions[:nlines-1]:
        data = mech('(%s)[%s]' % (rxn, species))
        options['color'] = chem_colors.pop()
        options['label'] = re.compile('\d+.\d+\*').sub('', str(mech(rxn).sum())).replace(' ', '')
        plot_date(date_objs, data.repeat(2,0), **options)

    data = zeros(data.shape, dtype = data.dtype)
    for rxn in reactions[nlines-1:]:
        data += mech('(%s)[%s]' % (rxn, species))
    options['color'] = chem_colors.pop()
    options['label'] = 'other'
    plot_date(date_objs, data.repeat(2,0), **options)
    
    options['color'] = 'black'
    options['marker'] = 'o'
    options['label'] = 'Chem'
    try:
        data = mech('%s[%s]'  % (conf.chem, species)).array()
    except:
        warn('Using sum of reactions for %(species)s' % locals())
        data = mech.make_net_rxn(species_obj, species_obj, False)[species_obj]

    plot_date(date_objs, data.repeat(2,0), **options)

    xlabel('Time')
    ylabel(units)
    ax.xaxis.set_major_formatter(DateFormatter('%H'))
    fig.autofmt_xdate()
    
    legend(loc=(0,-0.8), prop = FontProperties(size=10))
    return fig

def chem_plots(conf, nlines = 8, combine = [()], fmt = 'pdf'):
    mech = get_pure_mech('_'.join([conf.mechanism.lower(),conf.model.lower()]))
    tz = {'cmaq': 'UTC', 'camx': 'LST'}[conf.model.lower()]
    
    mrg_file = NetCDFFile(conf.mrgfile,'r')
    mech.set_mrg(mrg_file)
    
    date_objs = get_date(mrg_file, conf)

    for species, species_options in conf.species.iteritems():
        if species == 'NTR1':
            import pdb; pdb.set_trace()
        fig = chem_plot(conf, mech, date_objs, species, species_options, nlines = nlines, combine = combine)
        savefig(os.path.join(conf.outdir, '%s_IRR.%s' % (species, fmt)), format = fmt)

def phy_plot(conf, mech, date_objs, species, species_options):
    """
    conf - configuration obect that has title, 
           init, final, process, and species
           * title - string template that has access
                     to all local variables by name
           * mech - net_balance.Mechanism.Mechanism object
           * species - sp
    """
    units = mech.ipr.units
    fig = figure()
    grid(True)
    title(conf.title % locals())
    options = {}
    options.setdefault('c', 'k')
    options.setdefault('ls', '-')
    options.setdefault('lw', 3)
    options.setdefault('marker', 'o')
    data = mech('%s[%s]'  % (conf.init, species)).array()
    plot_date(date_objs[::2], data, **options)
    options.setdefault('label', 'Conc')
    options['marker'] = 'x'
    data = mech('%s[%s]'  % (conf.final, species)).array()
    plot_date(date_objs[1::2], data, **options)
    
    for process in conf.process.keys():
        options = conf.process[process]
        options.setdefault('linestyle', '-')
        options.setdefault('linewidth', 3)
        options.setdefault('label', process)
        options.setdefault('marker', 'None')
        options.setdefault('color', 'None')
        data = mech('(%s)[%s]' % (process, species)).array().repeat(2,0)
        if data.nonzero()[0].any() or not filter:
            plot_date(date_objs, data, **options)
    axis(**species_options)
    ax = axes()
    xlabel('Time')
    ylabel(units)
    ax.xaxis.set_major_formatter(DateFormatter('%H'))
    fig.autofmt_xdate()
    legend()
    return fig

def phy_plots(conf, filter = True, fmt = 'pdf'):
    mech = get_pure_mech('_'.join([conf.mechanism,conf.model]))
    mrg_file = NetCDFFile(conf.mrgfile,'r')
    mech.set_mrg(mrg_file)
    
    units = mech.ipr.units
    date_objs = get_date(mrg_file, conf)
    for species, species_options in conf.species.iteritems():
        try:
            fig = phy_plot(conf, mech,  date_objs, species, species_options)
            fig.savefig(os.path.join(conf.outdir, '%s_IPR.%s' % (species, fmt)), format = fmt)
        except KeyError, detail:
            warn(detail)
        
if __name__ == '__main__':
    conf = AttrDict()
    conf.mechanism = 'cb05'
    conf.model = 'camx'
    conf.mrgfile = '/Users/barronh/Development/net_balance/src/net_balance/testdata/test.mrg.nc'
    conf.species  = dict(HNO3=dict(ymin=None),
                        NO2=dict(ymin=None),
                        NO=dict(ymin=None),
                        NOx=dict(ymin=None),
                        NOyN=dict(ymin=None),
                        NTR=dict(ymin=None)
                        )
    conf.process = {'H_Trans': dict(color = 'b'), 
                    'V_Trans': dict(color = 'g'),
                    'Emissions': dict(color = 'r'),
                    'Deposit': dict(color = 'c'),
                    'Aero_Chem': dict(color = 'm'),
                    'CHEM': dict(color = 'y'),
                    'Motion': dict(color = 'teal'),
                    'TEMPADJ': dict(color = 'aqua')
                    }
    conf.init = 'INIT'
    conf.final = 'FCONC'
    conf.chem = 'CHEM'
    conf.title = '%(species)s IPR plot'
    conf.outdir = '.'
    conf.end_date = True
    phy_plots(conf)
#    chem_plot(conf, combine = [('RXN_01', 'RXN_02', 'RXN_03')])
    conf.title = '%(species)s IRR plot'
    chem_plots(conf)