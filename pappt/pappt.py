# !/usr/bin/env python -i
__doc__="""
PAPPT is the Process Analysis Post Processing Tools
These tools are composed of a driving function ext_mrg and
a helper class extracted.

ext_mrg fills an extracted object with rxn and process (for 
specified reactions and species) information according to the
shape, unit conversion, contribution, and normalizer variables 
specified.  The extracted object contains extracted information
that has been assigned to one of 5 boxes (NOCHG, VENT, VDET, 
HENT, HDET).

The extracted class provides an easy interface for categorizing
spatial values into their contribution categories. (see boxes 
above)  These categories can then be combined and normalized per
user specified values to retrieve the contribution to each process
"""

__all__ = ['MergedWriter', 'box_id', 'boxes', 'ext_mrg', 'extracted']


HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

import unittest
import operator, sys
from time import time

from numpy import vstack, ones, where, zeros, \
                  indices, squeeze, array, nansum, \
                  nan,nanmin, nanmax,hstack,newaxis, \
                  arange,dtype,isnan,ndarray
from pyPA.utils.sci_var import PseudoNetCDFFile,Pseudo2NetCDF
from pyPA.utils.timetuple import timeadd, timerange
from lagrangian import box_id,boxes

from pyPA.netcdf import NetCDFFile as ncf
    
def ext_mrg(pa_file,spc_iter,prc_iter,rxn_iter,shape=None,ipr_unitconversion=1,irr_unitconversion=1,ipr_contribution=1,irr_contribution=1,normalizer=1,kaxis=1):
    """
    ext_mrg is the drive horse and calls all the other functions and
    classes

    pa_file - must have dictionary property variables that returns
                something that can be converted to an array with 
                dimensions (time,layer,row,col)

    spc_iter, prc_iter, and rxn_iter - should iterate keys for the variables
                dictionary for species, processes, and reactions

    shape - 4 dimensional array(time,layer,row,col) with 0 and 1 values where 1s 
                define the shape of the analysis volume

    unit_conversion - supplies a factor to change units.  The factor must be broadcastable 
                to ipr and irr array dimensions

    contribution - supplies a factor to be applied to the cell that represents
                its contribution to the analysis volume

    normalizer - supplies an array that can be summed to be a normalization
                denominator

    \sigma_{ijk}{val_{ijk}*contribution_{ijk}*ucnv_{ijk}} /
    \sigma_{ijk}{normalization_{ijk}}
    """
    # timing the process for posterity
    startt=time()

    # If the shape is not defined, initialize it
    if shape==None:
        shape=ones(unit_conversion.shape,'f')
    
    # Improve efficiency by using envelope
    idx=indices(shape[0,:,:,:].shape).astype('f')+1
    idx*=shape[newaxis].max(1)
    idx[idx==0]=nan
    envelope=[
        slice(None),
        slice(nanmin(idx[0]-1).astype('i'),nanmax(idx[0]).astype('i')),
        slice(nanmin(idx[1]-1).astype('i'),nanmax(idx[1]).astype('i')),
        slice(nanmin(idx[2]-1).astype('i'),nanmax(idx[2]).astype('i'))
        ]
    shape=shape[envelope]
    try:
        normalizer=normalizer[envelope]
    except TypeError:
        if isinstance(normalizer,(int,float)):
            normalizer=shape[1:]*normalizer
        else:
            raise
        
    if type(irr_contribution) not in [float,int]:
        irr_contribution=irr_contribution[envelope]
    if type(irr_unitconversion) not in [float,int]:
        irr_unitconversion=irr_unitconversion[envelope]

    if type(ipr_contribution) not in [float,int]:
        ipr_contribution=ipr_contribution[envelope]
    if type(ipr_unitconversion) not in [float,int]:
        ipr_unitconversion=ipr_unitconversion[envelope]
    
    # Use shapes to define regions associated with (v and h) entrain,
    #  (v and h) detrain, initial, and hourly  data
    bxs=boxes(shape,kaxis)

    # create an extracted object to aggregate
    # moles
    agg=extracted(spc_iter,prc_iter,rxn_iter,bxs,normalizer)

    # Each reaction will be cut into appropriate boxes and summed to a single value
    # for each time.
    for ri,rxn_name in enumerate(rxn_iter):
        print >>sys.stderr, rxn_name
        rxn=array(pa_file.variables[rxn_name])[envelope]*irr_contribution*irr_unitconversion
        agg.aggregate_reaction(ri,rxn)
        
    # Each species/process combination will be cut into appropriate 
    # boxes and summed to a single value for each time.
    for si,spc_name in enumerate(spc_iter):
        # echo to user
        print >>sys.stderr, spc_name
        for pi,prc_name in enumerate(prc_iter):
            # for each variable
            var_name='_'.join((prc_name,spc_name))

            # echo to user
            # print >>sys.stderr, var_name,(si,pi)

            if var_name not in pa_file.variables.keys():
                print >>sys.stderr, "File does not contain %s it is assumed 0" % var_name
            else:    
                # get the variable value
                spc_prc=array(pa_file.variables[var_name])[envelope]*ipr_contribution*ipr_unitconversion

                # add it to the "stack"
                agg.aggregate_process(si,pi,spc_prc)

                # not strickly necessary, but using lots of memory
                try:
                    del spc_prc
                    del ipr.variables[var_name]
                except:
                    pass

    # print 'Time (min): ', (time()-startt)/60
    return agg

class extracted(object):
    """
    Extracted accumulates mol based process informaiton into 
    the separate boxes (horiz. (vertical) en(de)train, dilution)
    and has the functionality to convert the results to mrged
    and convert to ppm/ppb
    """
  
    def __init__(self,spc,prc,rxn,bxs,normalizer):
        """
        extracted takes dimensioning variables spc,prc,ntime
        and box definitions  (bxs), the air species (mols) and the
        keyword to id the initial process
        """
        # internalize boxes
        self.bxs=bxs
        self.spc=spc
        self.prc=list(prc)
        self.rxn=list(rxn)
        boxes=[(v,k) for k,v in box_id.iteritems()]; boxes.sort()
        boxes.pop(0)
        self.pseudo_proc=array(boxes)[:,1].tolist()
        self.pseudo_proc.append('EDHDIL')
        self.pseudo_proc.append('EDVDIL')
        self.pseudo_proc.append('TEMPADJ')
        self.prc.extend(self.pseudo_proc)
        ntime=bxs.shape[0]
        nboxes=bxs.shape[-1]
        # To ease computation, the normal shape of variables
        # will be coerced to add a 1 axis corresponding to boxes
        boxaxis=-1
        self.__new_shape=list(bxs.shape)
        self.__new_shape[boxaxis]=1

        # make it harder to edit
        self.__new_shape=tuple(self.__new_shape)

        # The air species is special since it will be sumed
        # separately.  The new shape will have single dimensioned
        # axes corresponding to spc and prc
        self.normalizer=normalizer.reshape(self.__new_shape)
        norm_shape=(ntime,nboxes)
        # norm_shape=list(bxs.shape)

        # The air species is summed within its boxes to give total box moles
        self.normalizer=nansum(nansum(nansum(self.normalizer*self.bxs,axis=-2),axis=-2),axis=-2).reshape(norm_shape)

        # Initialize data storage array
        self.process=zeros((ntime,len(spc),len(self.prc),nboxes),'f')
        self.reaction=zeros((ntime,len(rxn),nboxes),'f')
    
    def merge(self,prc={'INIT':[box_id.NOCHG,box_id.VDET,box_id.HDET]},convert=True,init_prc='INIT',final_prc='FCONC',exclude=['TEMPADJ','AVOL','UCNV']):
        """
        Prc is a dictionary of processes that maps to a list of 
        boxes to include in the merge process for that proc.  For
        processes not in that dictionary, entrained and unchanged
        will be used
        """

        # Default box set
        def_boxes=[box_id.NOCHG,box_id.VENT,box_id.HENT]

        # Dimension and create an output array
        ipr_output=zeros(self.process.shape[:-1],self.process.dtype)
        irr_output=zeros(self.reaction.shape[:-1],self.reaction.dtype)
        
        # For each process use the default box set if one is not
        # provided
        for pname in [p for p in self.prc if p not in self.pseudo_proc]:
            # get ordinal corresponding to dimension
            p=self.prc.index(pname)

            # get boxes to sum across
            boxes=prc.get(pname,def_boxes)
            
            # set values
            ipr_output[:,:,p]=nansum(self.process[:,:,p,boxes],axis=-1)/nansum(self.normalizer[:,boxes],axis=-1)[:,newaxis]
        
        # For pseudo-processes en(de)trainment, dilution and tempadj, over-write values
        # with delta in edtrain variable (INIT)
        init_id=self.prc.index(init_prc)
        
        #  Assumed order
        #  Pseudo Process: (1) V Detrain (2) V Entrain and Dilution (3) H Detrain (4) H Entrain and Dilution
        # 
        #  Equation subscripts u: nochg; vd: vertical detrainment; ve: vertical entrain and dilution; hd: horizontal detrainment; he: horizontal entrain and dilution
        #  1) (i_u+i_hd+i_vd)/(a_u+a_hd+a_vd)->(i_u+i_hd)/(a_u+a_hd)
        #  2) (i_u+i_hd)/(a_u+a_hd)->(i_u+i_hd+i_ve)/(a_u+a_hd+a_ve)
        #  3) (i_u+i_hd+i_ve)/(a_u+a_hd+a_ve)->(i_u+i_ve)/(a_u+a_ve)
        #  4) (i_u+i_ve)/(a_u+a_ve)->(i_u+i_ve+i_he)/(a_u+a_ve+a_he)
        i=lambda box: where(isnan(self.process),0,self.process)[:,:,init_id,box]
        a=lambda box: where(isnan(self.normalizer),0,self.normalizer)[:,box,newaxis]

        vdetrain=(
                    a(box_id.VDET)*(i(box_id.NOCHG)+i(box_id.HDET))-
                    i(box_id.VDET)*(a(box_id.NOCHG)+a(box_id.HDET))
                )/(
                    (a(box_id.NOCHG)+a(box_id.HDET))*(a(box_id.NOCHG)+a(box_id.HDET)+a(box_id.VDET))
                )

        ventrain=(
                    i(box_id.VENT)
                )/(
                    a(box_id.NOCHG)+a(box_id.HDET)+a(box_id.VENT)
                )
        
        vdilution=-(
                    (i(box_id.NOCHG)+i(box_id.HDET))*a(box_id.VENT)
                )/(
                    (a(box_id.NOCHG)+a(box_id.HDET))*(a(box_id.NOCHG)+a(box_id.HDET)+a(box_id.VENT))
                )

        hdetrain=(
                    a(box_id.HDET)*(i(box_id.NOCHG)+i(box_id.VENT))-
                    i(box_id.HDET)*(a(box_id.NOCHG)+a(box_id.VENT))
                )/(
                    (a(box_id.NOCHG)+a(box_id.VENT))*(a(box_id.NOCHG)+a(box_id.VENT)+a(box_id.HDET))
                )


        hentrain=(
                    i(box_id.HENT)
                )/(
                    a(box_id.NOCHG)+a(box_id.HENT)+a(box_id.VENT)
                )

        hdilution=-(
                    a(box_id.HENT)*(i(box_id.NOCHG)+i(box_id.VENT))
                )/(
                    (a(box_id.NOCHG)+a(box_id.VENT))*(a(box_id.NOCHG)+a(box_id.HENT)+a(box_id.VENT))
                )

        # Remove NANs for simplicity
        hdetrain=where(isnan(hdetrain),0,hdetrain)
        vdetrain=where(isnan(vdetrain),0,vdetrain)
        ventrain=where(isnan(ventrain),0,ventrain)
        hentrain=where(isnan(hentrain),0,hentrain)

        del i,a
        
        dilution=vdilution+hdilution
        
        # Assign values
        for pname,pid in [(k,v) for k,v in box_id.iteritems() if v!=0]:
            # get ordinal corresponding to dimension
            p=self.prc.index(pname)
            
            # get boxes to sum across
            boxes=prc.get(pname,def_boxes)

            # set values
            ipr_output[:,:,p]={'HENT':hentrain,'HDET':hdetrain, 'VDET': vdetrain,'VENT':ventrain}[pname]
        
        # Add dilution process
        p=self.prc.index('EDHDIL')
        ipr_output[:,:,p]=hdilution
        p=self.prc.index('EDVDIL')
        ipr_output[:,:,p]=vdilution
        
        #  Temperature adjustment accounts for difference between the final and sum of initial and process concentration
        #  TempAdj=Final-Init-dC
        procs=[self.prc.index(i) for i in self.prc if i not in [final_prc]+exclude]
        try:
            final_id=self.prc.index(final_prc)
            tempadj=ipr_output[:-1,:,final_id]-ipr_output[1:,:,init_id]
        except:
            tempadj=zeros(ipr_output.shape[:2],'f')[1:]

        # Add temperature adjustment process
        p=self.prc.index('TEMPADJ')
        ipr_output[1:,:,p]=tempadj

        irr_output[:,:]=nansum(self.reaction[:,:,def_boxes],-1)/nansum(self.normalizer[:,def_boxes],axis=-1)[:,newaxis]
        
        output=PseudoNetCDFFile()
        output.createDimension('TSTEP',ipr_output.shape[0])
        output.createDimension('SPECIES',ipr_output.shape[1])
        output.createDimension('PROCESS',ipr_output.shape[2])
        output.createDimension('RXN',irr_output.shape[1])
        
        v=output.createVariable('IPR','f',('TSTEP','SPECIES','PROCESS'))
        v.units="ppmV/time"
        v.assignValue(ipr_output)
        v=output.createVariable('IRR','f',('TSTEP','RXN'))
        v.units="ppmV/time"
        v.assignValue(irr_output)
        output.Process=''.join([i.ljust(16) for i in self.prc])
        output.Species=''.join([i.ljust(16) for i in self.spc])
        output.Reactions=''.join([i.ljust(16) for i in self.rxn])
        return output

    def boxagg(self,vals):
        cut=(array(vals).reshape(self.__new_shape)*self.bxs)
        return nansum(nansum(nansum(cut,axis=-2),axis=-2),axis=-2)
    
    def aggregate_reaction(self,rxn,vals):
        #  vals dim(time,layer,col,row)
        #  bxs dim(time,layer,col,row,bxs)
        self.reaction[:,rxn,:]=self.boxagg(vals)
    
    def aggregate_process(self,spc,prc,vals):
        #  vals dim(time,layer,col,row)
        #  bxs dim(time,layer,col,row,bxs)
        self.process[:,spc,prc,:]=self.boxagg(vals)

def MergedWriter(outpath,ipr_irr,shape,tflag):
    """
    Persists 3 arrays and meta data
    irr - 2 dimensional array (TSTEP, Reaction)
    ipr - 3 dimensional array (TSTEP, Species, Process)
    shape - 4  dimensional array (TSTEP,LAY,ROW,COL)
    TFLAG - 3 dimensional array (TSTEP,VAR,DATE-TIME)
    """
    ipr_irr.createDimension('VAR',3)
    ipr_irr.createDimension('DATE-TIME',2)
    shape_out=ipr_irr.createVariable('SHAPE','i',('TSTEP_STAG','LAY','ROW','COL'))
    shape_out.assignValue(shape)
    shape_out.units='onoff'
    shape_out.var_desc='SHAPE'.ljust(16)
    shape_out.long_name='SHAPE'.ljust(16)
    time=ipr_irr.createVariable('TFLAG','i',('TSTEP','VAR','DATE-TIME'))
    time[:,:,:]=tflag[:,newaxis,:]
    time.units='YYYYJJJ,HHDDMM'
    time.var_desc='TFLAG'.ljust(16)
    time.long_name='TFLAG'.ljust(16)
    
    
    return Pseudo2NetCDF().convert(ipr_irr,outpath)
    
class TestExtractor(unittest.TestCase):
    def setUp(self):
        spc=['NO','NO2','O3']
        prc=['INIT','ADV','DIF','EMIS','CHEM','FCONC']
        rxn=['IRR_1','IRR_2','IRR_3','IRR_4']
        ntime=4
        shape=zeros((5,3,4,5),'f')
        shape[0,0,0,0]=1
        shape[1,:2,0,0]=1
        shape[2,:2,:2,0]=1
        shape[3,:2,:2,:2]=1
        shape[4,:2,:2,:2]=1
        bxs=boxes(shape)
        self.shape=shape
        self.ext=extracted(spc,prc,rxn,bxs,ones(bxs.shape[:-1],'f'))
        self.vals=arange(240,dtype='f').reshape(4,3,4,5)
        self.ans=array([[0.,20.,0.,0.,0.],[140.,0.,0.,150.,0.],[530.,0.,0.,534.,0.],[1544.,0.,0.,0.,0.]],dtype='f')

    def testBoxAgg(self):
        self.assert_((self.ext.boxagg(self.vals)==self.ans).all())

    def testAggRxn(self):
        rxn=0
        self.ext.aggregate_reaction(rxn,self.vals)
        self.assert_((self.ext.reaction[:,rxn,:]==self.ans).all())

    def testAggPrc(self):
        spc=0;prc=1
        self.ext.aggregate_process(spc,prc,self.vals)
        self.assert_((self.ext.process[:,spc,prc,:]==self.ans).all())
    
    def testMerge(self):
        from lagrangian import box_id
        s=zeros((3,2,2,2),'i') # 3 times; 2x2x2 spatial grid
        s[0,0,0,0]=1 # turn on time 0, 0,0,0 cell
        s[1,...]=1 # turn on time 1 all cells
        s[2,1,1,1]=1 #  turn on time 2 1,1,1 cell
        
        v=zeros((2,2,2,2),'f') # create values for testing
        v[0,...]=arange(1,9,dtype='f').reshape((2,2,2)) # assign array with unique values
        v[1,...]=arange(1,9,dtype='f').reshape((2,2,2))*2
        
        # Create an extracted object with 3 procs, 1 species, 1 reaction, 2 times, boxes, and an
        # even normalizer
        ext=extracted(['O3'],['INIT','OTHER','FCONC'],['IRR'],boxes(s),ones((2,2,2,2),'f'))
        
        # Add values to reactions and both processes
        ext.aggregate_reaction(0,v)
        ext.aggregate_process(0,0,v)
        ext.aggregate_process(0,1,v)
        ext.aggregate_process(0,2,v*2)
        
        # Merge values
        mrg=ext.merge()
        # make irr and ipr values easily accessible
        irr=mrg.variables['IRR']
        ipr=mrg.variables['IPR']
        
        # Check time 0 values
        self.assertEqual(ipr[0,0,ext.prc.index('INIT')],1.) # INIT
        self.assertEqual(ipr[0,0,ext.prc.index('OTHER')],36./8.) # OTHER
        self.assertEqual(ipr[0,0,ext.prc.index('FCONC')],72./8.) # FCONC
        self.assertEqual(ipr[0,0,ext.prc.index('VENT')],5./2.) # VENT
        self.assertEqual(ipr[0,0,ext.prc.index('VDET')],0.) # VDET
        self.assertEqual(ipr[0,0,ext.prc.index('HENT')],30./8.) # HENT
        self.assertEqual(ipr[0,0,ext.prc.index('HDET')],0.) # HDET
        self.assertEqual(ipr[0,0,ext.prc.index('EDVDIL')],-0.5) # Dilution
        self.assertEqual(ipr[0,0,ext.prc.index('EDHDIL')],-2.25) # Dilution
        self.assertEqual(ipr[0,0,ext.prc.index('TEMPADJ')],0.) # Dilution

        # Check time 1 values
        self.assertEqual(ipr[1,0,ext.prc.index('INIT')],72./8.) # INIT
        self.assertEqual(ipr[1,0,ext.prc.index('OTHER')],16.) # other
        self.assertEqual(ipr[1,0,ext.prc.index('FCONC')],32.) # final
        self.assertEqual(ipr[1,0,ext.prc.index('VENT')],0) # VENT
        self.assertEqual(ipr[1,0,ext.prc.index('VDET')],array(0.14285715,'f')) # VDET
        self.assertEqual(ipr[1,0,ext.prc.index('HENT')],0) # HENT
        self.assertEqual(ipr[1,0,ext.prc.index('HDET')],array(6.85714293,'f')) # HDET
        self.assertEqual(ipr[1,0,ext.prc.index('EDVDIL')],0./8.) # Dilution
        self.assertEqual(ipr[1,0,ext.prc.index('EDHDIL')],0./8.) # Dilution
        self.assertEqual(ipr[1,0,ext.prc.index('TEMPADJ')],0./8.) # Dilution

        # Use large initial value to create dilution for testing
        v[0,0,0,0]=10.
        
        # recreate extract and remerge
        ext=extracted(['O3'],['INIT','OTHER','FCONC'],['IRR'],boxes(s),ones((2,2,2,2),'f'))
        ext.aggregate_reaction(0,v)
        ext.aggregate_process(0,0,v)
        ext.aggregate_process(0,1,v)
        ext.aggregate_process(0,2,v*2)
        mrg=ext.merge()
        irr=mrg.variables['IRR']
        ipr=mrg.variables['IPR']
        
        # Check values for time 1
        self.assertEqual(ipr[0,0,ext.prc.index('INIT')],10.) # INIT
        self.assertEqual(ipr[0,0,ext.prc.index('OTHER')],5.625) # OTHER
        self.assertEqual(ipr[0,0,ext.prc.index('FCONC')],11.25) # FCONC
        self.assertEqual(ipr[0,0,ext.prc.index('VENT')],2.5) # VENT
        self.assertEqual(ipr[0,0,ext.prc.index('VDET')],0) # VDET
        self.assertEqual(ipr[0,0,ext.prc.index('HENT')],3.75) # HENT
        self.assertEqual(ipr[0,0,ext.prc.index('HDET')],0) # HDET
        self.assertEqual(ipr[0,0,ext.prc.index('EDVDIL')],-5.) # Dilution
        self.assertEqual(ipr[0,0,ext.prc.index('EDHDIL')],-5.625) # Dilution
        self.assertEqual(ipr[0,0,ext.prc.index('TEMPADJ')],0./8.) # Dilution

    def runTest(self):
        pass


if __name__=='__main__':
    unittest.main()
