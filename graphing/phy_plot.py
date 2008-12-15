#!/usr/bin/env python
HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

#site-packages
import os,sys
from pynetcdf import NetCDFFile as ncdf
from numpy import array,newaxis,zeros,arange
from yaml import load
from datetime import datetime,timedelta
import time
from math import fmod
from pychartdir import *
from mx import DateTime
import calendar as cal
import pdb

__all__ = ['phy_plot']

# Object that takes a yaml file with several parameters which include time_slice, model, time_offsets, 
# species and process desired, as well as input data.
def phy_plot(job, ncfvariable= None):

    #Checks which model is used and and load default template
    phy_path = os.path.dirname(os.path.realpath(__file__))
    if job['model'].lower() == 'camx':
        def_y=load(file(os.path.join(phy_path,'phy_yaml/camx_template.yaml'),'r'))
    else:
        def_y=load(file(os.path.join(phy_path,'phy_yaml/cmaq_template.yaml'),'r'))
    
    # Check all user parameters include default paramaters if not user specified
    for k in def_y:
        job.setdefault(k,def_y[k])
    for spci,spc in enumerate(job['Species']):
        if len(spc)==2:
            job['Species'][spci].append([1]*len(spc[1]))
    
    inputpath = job['input_path'] 
    outputpath = job['output_path']
    
    episode_filename = os.path.basename(job['input_path'])
    
    ncf = ncdf(inputpath, 'r')
    list_species=[i.strip() for i in array(ncf.Species,ndmin=1).view('|S16')]
    prc_name_ncf = [i.strip() for i in array(ncf.Process,ndmin=1).view('|S16')]
    color=job['color']
    
    for spi,(spcn,spc_list,spc_wt) in enumerate(job['Species']):
            job['Species'][spi]=[spcn,[list_species.index(i) for i in spc_list],spc_wt]
    for pri,(prcn,prc_list) in enumerate(job['Process']):
        job['Process'][pri]=[prcn,[prc_name_ncf.index(i) for i in prc_list]]
    
    ipr = ncf.variables['IPR'] 
    timetest = ncf.variables['TFLAG'][:,0,:].tolist()

    for ni,i in enumerate(timetest):
        timetest[ni][-1]= str(timetest[ni][-1]/10000).zfill(2)
        
    sidx=timetest.index(job['time_slice'][0]);eidx=timetest.index(job['time_slice'][1])
    timetest = timetest[sidx:eidx+1]
    
    timetemp=[]
    for ni,i in enumerate(timetest):
        timetemp.append(DateTime.ISO.ParseDateTime(time.strftime("%Y-%m-%d %H:%M:%S", time.strptime(str(timetest[ni]), "[%Y%j, '%H']"))))
    
    
    dataX = range(sidx,eidx+1)
    if job['begin_hour'] == True:
        timetemp.append(timetemp[-1]+.041666666666666666)
    else:
        timetemp.insert(0,timetemp[0]-.041666666666666666)
    data_ins = range(len(dataX))
    data_step = range(len(dataX)+1)
    
    timehours1=[(timetemp[ni]+(.041666666666666666*job['display_time_offset'])).hour for ni,i in enumerate(timetemp)]
    tt=[timetemp[ni]+(.041666666666666666*job['display_time_offset']) for ni,i in enumerate(timetemp)]
    begin ='%s %02d' % (cal.month_abbr[tt[0].month], tt[0].day); end ='%s %02d' % (cal.month_abbr[tt[-1].month], tt[-1].day)
    label=['%02d' % i for i in timehours1]
    label[0] = label[0]+'\n'+begin
    label[-1] = label[-1]+'\n'+end
    
    DATE=ncf.SDATE[0]
    date=str(datetime(DATE/1000,1,1)+timedelta(days=fmod(DATE,1000.)-1)).split()
    DATE=date[0].replace('-','')

    for spi,(spcn,spc_list,spc_wt) in enumerate(job['Species']):
        # Create a XYChart object of size 500 x 270 pixels, with white
        # background, black border, 1 pixel 3D border effect and rounded corners
        c = XYChart(800, 450, 0xffffff, 0x000000, 1)
        c.setRoundedFrame()
        
        # Set the plotarea at (55, 60) and of size 520 x 200 pixels, with white (ffffff)
        # background. Set horizontal and vertical grid lines to grey (cccccc).
        c.setPlotArea(60, 90, 700, 290, 0xffffff, -1, -1, 0xcccccc, 0xffffff)
        
        # Add a legend box at (55, 32) (top of the chart) with horizontal layout. Use 9 pts
        # Arial Bold font. Set the background and border color to Transparent.
        legendBox = c.addLegend2(140, 32, 4, "arialbd.ttf", 11)
        legendBox.setBackground(Transparent, Transparent)
        
        # Add a title box to the chart using 15 pts Times Bold Italic font. The text is white
        # (ffffff) on a deep blue (000088) background, with soft lighting effect from the
        # right side.
        c.addTitle('%s %s' % (spcn, date[0]), "arial.ttf", 16, 0x000000).setBackground(0xffffff)
        
        # Set the labels on the x axis.
        c.xAxis().setLabels(label)
        
        # Add a title to the y axis
        c.yAxis().setTitle(job['units_Conc'])
        
        if job['Scale'] == True:
            c.yAxis().setLinearScale(job['Scale_lowerLimit'],job['Scale_upperLimit'])
        
        # Add a title to the x axis
        c.xAxis().setTitle(job['units_Time'])
        
        # Add title to the bottom of the chart using 7.5pt Arial font. The text is black 0x000000.
        c.addTitle2(Bottom, episode_filename, "arial.ttf", 10, 0x000000,0xffffff)   
    
        step_layer= range(len(job['Process']))
        prc_ini = str()
        for pri,(prcn,prc_list) in enumerate(job['Process']):
            spc_a=array(ipr)[dataX,:,:][:,spc_list,:][:,:,prc_list]
            spc_a=spc_a*array(spc_wt)[newaxis,:,newaxis]
            spc_sum=spc_a.sum(1).sum(1)
            if prcn == 'Conc_Initial':
                instant_layer = c.addLineLayer2()
                instant_layer.addDataSet(list(spc_sum), color[pri], 'Conc_Initial').setDataSymbol(DiamondSymbol, 11)
                instant_layer.setXData(data_ins)
            else:
                temp=arange(0,spc_sum.size+1)
                temp[-1]=temp[-2]
                step_layer[pri] = c.addStepLineLayer(list(spc_sum[temp]), color[pri], prcn)
                step_layer[pri].setLineWidth(5)
                step_layer[pri].setXData(data_step)
            
        c.makeChart('%s/phy_plot_%s_%s.png' %(outputpath,spcn,DATE))

if __name__ == '__main__':
    if sys.argv < 1:
        print '*****Error: Missing Argument****'
        print '**To Get Copy of yaml input template set Copytemplates as first argument**'
        print '***and to choose where to store templates -i outpath for templates***'
    #TEST CASE
    elif sys.argv[1] == 'Testcase': 
        from pyPA.utils.sci_var import PseudoNetCDFFile as pncf, PseudoNetCDFVariable as pncv
        from numpy import array
        testfile=pncf()
        testfile.createDimension('SPECIES',1)
        testfile.createDimension('TSTEP',24)
        testfile.createDimension('PROCESS',18) 
        test = testfile.createVariable('IPR','f',('TSTEP','SPECIES','PROCESS'))
        testfile.createDimension('VAR',3)
        testfile.createDimension('DATE-TIME',2)
        times = testfile.createVariable('TFLAG','i',('TSTEP','VAR','DATE-TIME'))
        testfile.Species="TESTCASE        "
        testfile.Process= 'INIT            CHEM            EMIS            XADV            YADV            ZADV            ADJC            HDIF            VDIF            DDEP            FCONC           VENT            VDET            HENT            HDET            EDHDIL          EDVDIL          TEMPADJ         '
        testfile.SDATE= [2006229]
        tmptime=[[[2006229,   10000],         [2006229,   10000],         [2006229,   10000]],         [[2006229,   20000],         [2006229,   20000],         [2006229,   20000]],         [[2006229,   30000],         [2006229,   30000],         [2006229,   30000]],         [[2006229,   40000],         [2006229,   40000],         [2006229,   40000]],         [[2006229,   50000],         [2006229,   50000],         [2006229,   50000]],         [[2006229,   60000],         [2006229,   60000],         [2006229,   60000]],         [[2006229,   70000],         [2006229,   70000],         [2006229,   70000]],         [[2006229,   80000],         [2006229,   80000],         [2006229,   80000]],         [[2006229,   90000],         [2006229,   90000],         [2006229,   90000]],         [[2006229,  100000],         [2006229,  100000],         [2006229,  100000]],         [[2006229,  110000],         [2006229,  110000],         [2006229,  110000]],         [[2006229,  120000],         [2006229,  120000],         [2006229,  120000]],         [[2006229,  130000],         [2006229,  130000],         [2006229,  130000]],         [[2006229,  140000],         [2006229,  140000],         [2006229,  140000]],         [[2006229,  150000],         [2006229,  150000],         [2006229,  150000]],         [[2006229,  160000],         [2006229,  160000],         [2006229,  160000]],         [[2006229,  170000],         [2006229,  170000],         [2006229,  170000]],         [[2006229,  180000],         [2006229,  180000],         [2006229,  180000]],         [[2006229,  190000],         [2006229,  190000],         [2006229,  190000]],         [[2006229,  200000],         [2006229,  200000],         [2006229,  200000]],         [[2006229,  210000],         [2006229,  210000],         [2006229,  210000]],         [[2006229,  220000],         [2006229,  220000],         [2006229,  220000]],         [[2006229,  230000],         [2006229,  230000],         [2006229,  230000]],         [[2006230,       0],         [2006230,       0],         [2006230,       0]]]
        for i in range(len(tmptime)):
         for p in range(len(tmptime[i][:][:])):
          for n,ni in enumerate(tmptime[i][p][:]):
           times[i][p][n] = float(ni)
        tmptest=[[[  4.15039092e-01,  -3.46339035e+00,   3.09677649e+00,           2.55673174e-02,  -1.81295410e-01,   4.31699343e-02,          -1.26131612e-03,  -3.99496173e-03,   8.34290013e-02,          -1.80426752e-03,   1.01428427e-01,   0.00000000e+00,           8.91927108e-02,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -0.00000000e+00,   2.23517418e-07]],        [[  1.01448745e-01,  -8.56803894e+00,   8.52479172e+00,          -3.14713153e-03,  -1.56799674e-01,   2.15585381e-02,          -2.28288327e-03,  -1.52531378e-02,   5.01981378e-02,          -2.96872854e-03,   1.09970465e-03,   0.00000000e+00,           5.15931286e-02,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -0.00000000e+00,  -3.32598574e-07]],        [[  1.10103469e-03,  -1.00336504e+01,   1.02019167e+01,          -3.88754085e-02,  -1.65796980e-01,   3.10196169e-02,          -2.97858729e-03,  -2.00538058e-02,   3.84413227e-02,          -3.21150618e-03,   8.16731621e-03,   0.00000000e+00,           2.56840984e-04,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -0.00000000e+00,  -1.97440386e-06]],        [[  8.17066431e-03,  -9.48297882e+00,   9.72537708e+00,          -8.25400129e-02,  -1.48804754e-01,   2.44857855e-02,          -3.01884767e-03,  -2.11105663e-02,   2.26942990e-02,          -2.90278182e-03,   3.96084934e-02,   0.00000000e+00,           2.36360152e-04,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -0.00000000e+00,  -3.35276127e-07]],        [[  3.96137238e-02,  -8.94899178e+00,   9.44221497e+00,          -1.00060716e-01,  -2.44391650e-01,   1.55646401e-02,          -2.72348407e-03,  -2.38108505e-02,  -1.00003732e-02,          -2.59389542e-03,   1.65428743e-01,   0.00000000e+00,           6.08021393e-04,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -0.00000000e+00,   1.19209290e-07]],        [[  1.65423319e-01,  -8.44221210e+00,   8.90409088e+00,          -9.32137221e-02,  -2.28767931e-01,   1.19806416e-02,          -2.33306014e-03,  -2.40151156e-02,   8.28571431e-03,          -2.51914375e-03,   2.63074458e-01,   1.29534442e-06,          -3.17847505e-02,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -1.86137552e-03,  -2.38418579e-07]],        [[  2.63088763e-01,  -8.33423805e+00,   9.25724983e+00,          -1.36853531e-01,  -4.91358519e-01,   1.62423197e-02,          -2.94986204e-03,  -3.62655781e-02,  -7.54895359e-02,          -3.01825674e-03,   5.27233601e-01,   7.44385794e-02,           1.22666208e-03,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -4.83972579e-03,   5.96046448e-08]],        [[  5.27299523e-01,  -8.06521130e+00,   9.07385921e+00,          -2.42319569e-01,  -6.70686185e-01,   1.47283180e-02,          -2.19775178e-03,  -4.18188758e-02,   6.30116910e-02,          -3.97608522e-03,   5.01196802e-01,   2.02787135e-07,          -1.49687856e-01,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -1.80497440e-03,   5.36441803e-07]],        [[  5.01231670e-01,  -8.02879143e+00,   9.30801773e+00,          -4.08534378e-01,  -6.44423604e-01,   8.82820785e-03,          -2.64116609e-03,  -6.82353154e-02,   8.77440423e-02,          -5.33439685e-03,   7.72850573e-01,   0.00000000e+00,           2.49884315e-02,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -0.00000000e+00,   7.74860382e-07]],        [[  7.72812605e-01,  -8.88064480e+00,   1.06944666e+01,          -5.16222656e-01,  -6.68560326e-01,   6.17920095e-03,          -3.02923494e-03,  -9.73219872e-02,   1.23567268e-01,          -8.93061329e-03,   1.47222853e+00,   0.00000000e+00,           4.99132499e-02,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -0.00000000e+00,  -1.07288361e-06]],        [[  1.47198665e+00,  -1.10186596e+01,   1.38120079e+01,          -6.21128619e-01,  -6.76825166e-01,   8.88649188e-03,          -1.30916457e-03,  -1.25116900e-01,   2.19126388e-01,          -1.78512130e-02,   3.06925321e+00,   1.01737129e-02,           1.60317682e-02,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -8.06835573e-03,  -1.43051147e-06]],        [[  3.06856823e+00,  -1.07044106e+01,   1.70428123e+01,          -7.52368033e-01,  -7.49084711e-01,   3.15314010e-02,           4.08147881e-03,  -1.77031919e-01,   3.80629689e-01,          -3.98479030e-02,   8.16418552e+00,   7.58502781e-02,           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -1.65453367e-02,   9.53674316e-07]],        [[  8.16080284e+00,  -8.41549015e+00,   1.53915939e+01,          -5.57414234e-01,  -8.64693284e-01,  -9.02693532e-03,           2.21665669e-03,  -2.05351561e-01,   1.87457752e+00,          -5.99057227e-02,   1.44638491e+01,   5.94519198e-01,           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -1.44798040e+00,   9.53674316e-07]],        [[  1.44566307e+01,  -6.17867947e+00,   8.32131577e+00,           1.70661777e-01,  -3.81479651e-01,   1.27887040e-01,          -2.40133125e-02,  -1.07111074e-01,   3.19343233e+00,          -3.59601825e-02,   1.43554296e+01,   1.34739792e+00,           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -6.53465271e+00,  -9.53674316e-07]],        [[  1.43493690e+01,  -3.54451132e+00,   2.06592536e+00,           1.30346954e-01,  -1.01087861e-01,   9.46561620e-02,          -8.67424812e-03,  -1.97707936e-02,   2.16237783e+00,          -3.94287053e-03,   3.52848864e+00,   4.40486938e-01,           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -1.20366859e+01,  -1.43051147e-06]],        [[  3.52699018e+00,  -1.51547396e+00,   1.88986969e+00,          -2.55647182e-01,  -5.55208586e-02,   3.43074016e-02,          -1.04896110e-02,  -1.44081796e-02,   3.03723961e-01,          -1.81535468e-03,   3.39542007e+00,   8.59080479e-02,           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -5.92024565e-01,   4.76837158e-07]],        [[  3.39434195e+00,  -2.09427476e+00,   1.69707012e+00,          -5.51970303e-01,  -1.18645914e-01,   2.11374499e-02,          -1.29356803e-02,  -1.28215160e-02,   2.38442138e-01,          -1.26803631e-03,   2.24950266e+00,   5.10454644e-03,           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -3.14677387e-01,   0.00000000e+00]],        [[  2.24849319e+00,  -2.32849956e+00,   1.45312095e+00,          -5.48735380e-01,  -1.68731585e-01,   1.56334769e-02,          -8.46210588e-03,  -5.18759154e-03,   3.33206505e-01,          -7.56656227e-04,   7.82909930e-01,   2.34926473e-02,           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -2.30664089e-01,   1.19209290e-07]],        [[  7.82546937e-01,  -1.83896542e+00,   1.37594986e+00,          -1.81181401e-01,  -1.21686630e-01,   1.32097788e-02,          -2.71477993e-03,  -1.93086115e-03,   3.36707115e-01,          -4.91592276e-04,   3.19240570e-01,   1.96875166e-03,           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -4.41708751e-02,  -3.57627869e-07]],        [[  3.19118351e-01,  -1.51102841e+00,   1.23436975e+00,          -4.71397080e-02,  -8.28588009e-02,   1.47549333e-02,          -1.39944284e-04,  -1.91248942e-03,   3.37479085e-01,          -4.42206394e-04,   2.41647273e-01,   2.92139454e-03,           0.00000000e+00,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -2.34744232e-02,  -2.98023224e-07]],        [[  2.41553724e-01,  -1.38673961e+00,   1.16866922e+00,          -4.11386304e-02,  -8.28774124e-02,   1.77893378e-02,           5.55822917e-04,  -2.02063145e-03,   3.36533427e-01,          -4.60448849e-04,   2.50401199e-01,   3.78527853e-04,           3.81796708e-04,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -2.22407863e-03,   2.38418579e-07]],        [[  2.50342846e-01,  -1.37092018e+00,   1.18228436e+00,          -4.50761020e-02,  -1.03670850e-01,   2.08732933e-02,           9.47726599e-04,  -2.34447513e-03,   3.47631872e-01,          -5.39905042e-04,   2.80932903e-01,   2.10231447e-04,           1.97517057e-03,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -7.81176670e-04,   8.94069672e-08]],        [[  2.80945688e-01,  -1.39174628e+00,   1.26298344e+00,          -5.16617522e-02,  -1.58258229e-01,   2.52211541e-02,          -2.36925425e-05,  -3.02077457e-03,   3.07481110e-01,          -6.66284119e-04,   2.90701270e-01,   0.00000000e+00,           1.94468666e-02,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -0.00000000e+00,   5.96046448e-08]],        [[  2.90820181e-01,  -1.65909207e+00,   1.48323274e+00,          -4.44809496e-02,  -1.36264279e-01,   2.08873693e-02,           3.12612276e-04,  -3.97682050e-03,   2.15136856e-01,          -8.94676370e-04,   1.98358759e-01,   0.00000000e+00,           3.26778553e-02,   0.00000000e+00,   0.00000000e+00,          -0.00000000e+00,  -0.00000000e+00,  -1.19209290e-07]]]
        for i in range(len(tmptest)):
         for p in range(len(tmptest[i][:][:])):
          for n,ni in enumerate(tmptest[i][p][:]):
           test[i][p][n] = float(ni)
        pdb.set_trace()
        phy_path = os.path.dirname(os.path.realpath(__file__))
        testcase_path = os.path.join(phy_path,'phy_yaml/Testcase.yaml')
        job=load(file(testcase_path,'r'))
        phy_plot(job, testfile)
    elif sys.argv[1] == 'Copytemplates':
        phy_path = os.path.dirname(os.path.realpath(__file__))
        yaml_path_camx=os.path.join(phy_path,'phy_yaml/camx_template.yaml')
        yaml_path_cmaq=os.path.join(phy_path,'phy_yaml/cmaq_template.yaml')
        from shutil import copy
        copy(yaml_path_camx,'./camx_template.yaml')
        copy(yaml_path_cmaq,'./cmaq_template.yaml')
    elif sys.argv < 1:
        print '*****Error: Missing Argument****'
        print '**To Get Copy of yaml input template set Copytemplates as first argument**'
        print '***and to choose where to store templates -i outpath for templates***'
    else:
        job=load(file(sys.argv[1],'r')) #argument file
        phy_plot(job)
