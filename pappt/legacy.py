HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = "$RevisionNum$"

"""
Legacy is a module for converting new pappt outputs to old text based formats
"""

try:
	from pynetcdf import NetCDFFile as ncf
except:
	from Scientific.IO.NetCDF import NetCDFFile as ncf
from numpy import array

def LegacyMerged(outpath,inpath):
    """
    This function takes a CAMx based pappt output file (netcdf format) and converts
    the data into a text file that net_balance can read and use
    """
    if type(inpath)==str:
        inf=ncf(inpath)
    else:
        inf=inpath
    if type(outpath)==str:
        out=file(outpath,'w')
    else:
        out=outpath
    print >> out, '"%s' % inf.iprfile
    print >> out, '"%s' % inf.irrfile
    spcs=[s.strip() for s in array(inf.Species,ndmin=1).view('S16')]
    prcs=[p.strip() for p in array(inf.Process,ndmin=1).view('S16')]
    for ti,(d,t) in enumerate(inf.variables['TFLAG'][:,0,:]):
        print >> out, "Time =%06d" % (t*(10000/array(inf.TSTEP)),)
        print >> out, '!"Rxn No"    "Int Rate"'
        rr=inf.variables['IRR'][ti,:]
        for ri,r in enumerate(rr):
            print >> out, '{%03s}    %13.5e' % (ri+1,r)
        print >> out, ";"
        print >> out, "! Species      Initial conc.    Chemistry        Area emi.        Pt source emi.   PiG change       West b. adv.     East b. adv.     South b. adv.    North b. adv.    Bottom b. adv.   Top b. adv.      Dil. in the vert West b. diff.    East b. diff.    South b. diff.   North b. diff.   Bottom b. diff.  Top b. diff.     Dry dep.         Wet dep.         Aerosol chemistr Dilution         En(De)trainment  Final conc.      Units conversion Average cell vol"
        order2match=['INIT','CHEM','EMIS','PTEMIS','PIG','A_W','A_E','A_S','A_N','A_B','A_T','DIL','D_W','D_E','D_S','D_N','D_B','D_T','DDEP','WDEP','AERCHEM','EDDIL','EDTRAIN','FCONC','UCNV','AVOL']
        pr=array(inf.variables['IPR'])[ti,:,:]
        for si,spc in enumerate(spcs):
            pvals=[]
            for p in order2match:
                if p=='EDTRAIN':
                    val=pr[si,[prcs.index('HENT'),prcs.index('HDET'),prcs.index('VDET'),prcs.index('VENT')]].sum()
                elif p=='AERCHEM':
                    try:
                        val=pr[si,prcs.index('INORGACHEM')]+pr[si,prcs.index('ORGACHEM')]+pr[si,prcs.index('AQACHEM')]
                    except:
                        try:
                            val=pr[si,prcs.index(p)]
                        except:
                            val=0
                elif p=='EDDIL':
                    val=pr[si,[prcs.index('EDHDIL'),prcs.index('EDVDIL')]].sum()
                else:
                    try:
                        val=pr[si,prcs.index(p)]
                    except:
                        val=0
                pvals.append(val)
            pvals=tuple(pvals)
            print >> out, ('"%-8s" '+'%14.6e   '*len(pvals)) % ((spc,)+pvals)
        print >> out, ";"
    print >> out, "|"
        
    
def LegacyMergedCMAQ(outpath,inpath):
    """
    This function takes a CMAQ based pappt output file (netcdf format) and converts
    the data into a text file that net_balance can read and use
    
    The most notable difference from LegacyMerged is that the reactions must be
    reordered to account for the differences in CAMx and CMAQ implementations of 
    CB4
    """
    if type(inpath)==str:
        inf=ncf(inpath)
    else:
        inf=inpath
    if type(outpath)==str:
        out=file(outpath,'w')
    else:
        out=outpath
    
    r2r={
      1: [1], 
      2: [2], 
      3: [3], 
      4: [4], 
      5: [5], 
      6: [6], 
      7: [7], 
      8: [8], 
      9: [9], 
      10: [10,11], 
      11: [12],
      12: [13],
      13: [14],
      14: [15],
      15: [16],
      16: [17],
      17: [18],
      18: [19],
      19: [20],
      20: [21],
      21: [22],
      22: [23],
      23: [24],
      24: [25],
      25: [26],
      26: [27],
      27: [28],
      28: [29],
      29: [30],
      30: [31],
      31: [32],
      32: [33],
      33: [34],
      34: [35],
      35: [36],
      36: [37],
      37: [38],
      38: [39],
      39: [40],
      40: [41],
      41: [42],
      42: [43],
      43: [44],
      44: [45],
      45: [46],
      46: [47],
      47: [48],
      48: [49],
      49: [50],
      50: [51],
      51: [52],
      52: [53],
      53: [54],
      54: [55],
      55: [56],
      56: [57],
      57: [58],
      58: [59],
      59: [60],
      60: [61],
      61: [62],
      62: [63],
      63: [64],
      64: [65],
      65: [66],
      66: [67],
      67: [68],
      68: [69],
      69: [72],
      70: [71],
      71: [73],
      72: [70],
      73: [74],
      74: [75],
      75: [76],
      76: [77],
      77: [78],
      78: [79],
      79: [80],
      80: [81],
      81: [82],
      82: [83],
      83: [84],
      84: None,
      85: None,
      86: [85],
      87: [86],
      88: [87],
      89: [88],
      90: None,
      91: None,
      92: [89],
      93: [90],
      94: [91],
      95: [92],
      96: [93]
      }
    print >> out, '"%s' % inf.iprfile
    print >> out, '"%s' % inf.irrfile
    spcs=[s.strip() for s in array(inf.Species,ndmin=1).view('S16')]
    prcs=[p.strip() for p in array(inf.Process,ndmin=1).view('S16')]
    for ti,(d,t) in enumerate(inf.variables['TFLAG'][:,0,:]):
        print >> out, "Time =%06d" % (t*(10000/array(inf.TSTEP)),)
        print >> out, '!"Rxn No"    "Int Rate"'
        rr=array(inf.variables['IRR'][ti,:])
        for camxri in range(1,97):
            cmaqri=r2r[camxri]
            if cmaqri==None:
                rv=0.
            else:
                cmaqri=[id-1 for id in cmaqri]
                rv=rr[cmaqri].sum()
                
            print >> out, '{%03s}    %13.5e' % (camxri,rv)
            
        print >> out, ";"
        print >> out, "! Species      Initial conc.    Chemistry        Area emi.        Pt source emi.   PiG change       West b. adv.     East b. adv.     South b. adv.    North b. adv.    Bottom b. adv.   Top b. adv.      Dil. in the vert West b. diff.    East b. diff.    South b. diff.   North b. diff.   Bottom b. diff.  Top b. diff.     Dry dep.         Wet dep.         Aerosol chemistr Dilution         En(De)trainment  Final conc.      Units conversion Average cell vol"
        order2match=['INIT','CHEM','EMIS','PTEMIS','PIG','A_W','A_E','A_S','A_N','A_B','A_T','DIL','D_W','D_E','D_S','D_N','D_B','D_T','DDEP','WDEP','AERCHEM','EDDIL','EDTRAIN','FCONC','UCNV','AVOL']
        pr=array(inf.variables['IPR'])[ti,:,:]
        for si,spc in enumerate(spcs):
            pvals=[]
            for p in order2match:
                if p=='EDTRAIN':
                    val=pr[si,[prcs.index('HENT'),prcs.index('HDET'),prcs.index('VDET'),prcs.index('VENT')]].sum()
                elif p=='AERCHEM':
                    try:
                        val=pr[si,prcs.index('INORGACHEM')]+pr[si,prcs.index('ORGACHEM')]+pr[si,prcs.index('AQACHEM')]
                    except:
                        try:
                            val=pr[si,prcs.index(p)]
                        except:
                            val=0
                elif p=='EDDIL':
                    val=pr[si,[prcs.index('EDHDIL'),prcs.index('EDVDIL')]].sum()
                elif p=='A_W':
                    val=pr[si,prcs.index('XADV')]
                elif p=='A_S':
                    val=pr[si,prcs.index('YADV')]
                elif p=='A_T':
                    val=pr[si,prcs.index('ZADV')]
                elif p=='D_W':
                    val=pr[si,prcs.index('HDIF')]
                elif p=='D_T':
                    val=pr[si,prcs.index('VDIF')]
                elif p=='DIL':
                    val=pr[si,prcs.index('ADJC')]
                else:
                    try:
                        val=pr[si,prcs.index(p)]
                    except:
                        val=0
                pvals.append(val)
            pvals=tuple(pvals)
            print >> out, ('"%-8s" '+'%14.6e   '*len(pvals)) % ((spc,)+pvals)
        print >> out, ";"
    print >> out, "|"
