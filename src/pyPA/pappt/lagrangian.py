__all__ = ['box_id', 'boxes']

HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

import unittest
from numpy import indices,zeros,nanmax,nanmin,nansum,where,vstack,squeeze,nan
from ..utils.util import AttrDict

box_id=AttrDict({'NOCHG': 0, 'VENT': 1, 'VDET': 2, 'HENT': 3, 'HDET': 4})
def kvals(shape,kaxis):
    #get indices as an easy way to id k layer
    i=indices(shape.shape)
    #make one indexed for comparison to 0 cells
    i[kaxis,:,:,:,:]+=1

    #return values set to zero outside the shape
    return (i[kaxis,:,:,:,:])*where(shape==0,nan,shape)
  
def boxes(shape,kaxis=1):
    """boxes uses the now and next occupied cell matrices
    to determine which cells will be used for en(de)trainment, dillution
    and other processes
    """
    shape_now=shape[:-1,:,:,:]
    shape_next=shape[1:,:,:,:]

    #For debugging output
    debug=False

    #Initialize output array
    output=zeros(shape_now.shape+(len(box_id),),'bool_')

    #For easy math, all arrays will be dimensioned (time,layer,col,row)
    #column level information will be dimensioned (time,1,col,row)
    statshape=list(shape_now.shape)
    statshape[kaxis]=1

    #Get the old values
    okvals=kvals(shape_now,kaxis)
    okmax=nanmax(okvals,axis=kaxis).reshape(statshape)
    okmin=nanmin(okvals,axis=kaxis).reshape(statshape)
    ocount=nansum(shape_now,axis=kaxis).reshape(statshape)

    #Get the new
    nkvals=kvals(shape_next,kaxis)
    nkmax=nanmax(nkvals,axis=kaxis).reshape(statshape)
    nkmin=nanmin(nkvals,axis=kaxis).reshape(statshape)
    ncount=nansum(shape_next,axis=kaxis).reshape(statshape)

    #The instances where the new cell is higher than the old top
    highertop=where(nkvals>okmax,1,0)
    #The instances where the new cell is lower than the old bottom
    lowerbottom=where(nkvals<okmin,1,0)

    #The instances where the old cell is higher than the new top
    lowertop=where(okvals>nkmax,1,0)
    #The instances where the old cell is lower than the new bottom
    higherbottom=where(okvals<nkmin,1,0)

    #Instances where the column is unchanged
    same=(where(okmax==nkmax,1,0)+where(okmin==nkmin,1,0)+where(ocount==ncount,1,0))==3
    #The inverse of same
    notsame=(same==False).view('int8')

    #Column usage
    cololdyes=where(shape_now.sum(axis=kaxis).reshape(statshape)>0,1,0)
    colnewyes=where(shape_next.sum(axis=kaxis).reshape(statshape)>0,1,0)
    cololdno=(cololdyes==0).view('int8')
    colnewno=(colnewyes==0).view('int8')
    colboth=((where(cololdyes,1,0)+where(colnewyes,1,0))==2).view('int8')

    #Vertical entrainment is defined by the growth of an already
    #occupied column
    output[:,:,:,:,box_id.VENT]= \
        where((notsame+highertop+colboth)==3,True,False) + \
        where((notsame+colboth+lowerbottom)==3,True,False)
    if debug:
        print 'Hour ID\n', '%6i'*25 % tuple(range(25))
        print 'Count\n', '%6i'*25 % tuple(shape_next.sum(-1).sum(-1).sum(-1).tolist())
        print 've\n', '%6i'*25 % tuple(output[:,:,:,:,box_id.VENT].sum(-1).sum(-1).sum(-1).tolist())

    #Vertical detrainment is defined by the shrinkage of an already
    #occupied column
    output[:,:,:,:,box_id.VDET]= \
      where((lowertop+colboth+notsame)==3,True,False) + \
      where((higherbottom+colboth+notsame)==3,True,False)
    if debug:
        print 'vd\n', '%6i'*25 % tuple(output[:,:,:,:,box_id.VDET].sum(-1).sum(-1).sum(-1).tolist())
    #Horizontal Entrainment is defined by the use of a new column
    output[:,:,:,:,box_id.HENT]=where((cololdno+colnewyes+shape_next)==3,True,False)
    if debug:
        print 'he\n', '%6i'*25 % tuple(output[:,:,:,:,box_id.HENT].sum(-1).sum(-1).sum(-1).tolist())
  
    #Horizontal Detrainment is defined by no longer using a column
    output[:,:,:,:,box_id.HDET]=where((cololdyes+colnewno+shape_now)==3,True,False)
    if debug:
        print 'hd\n', '%6i'*25 % tuple(output[:,:,:,:,box_id.HDET].sum(-1).sum(-1).sum(-1).tolist())

    #No change: all other cells that are in use
    output[:,:,:,:,box_id.NOCHG]=(output.sum(-1)==0)*shape_now
    if debug:
        print 'No Change\n', '%6i'*25 % tuple(output[:,:,:,:,box_id.NOCHG].sum(-1).sum(-1).sum(-1).tolist())
        t=11;k=2;i=4;j=4
        key=(t,k,i,j)
        colkey=(t,0,i,j)
        print shape_now.shape
        print shape_now[key], shape_next[key]
        print 'NS',notsame[colkey],'CB',colboth[colkey]
        print 'CO', cololdyes[colkey], 'CN', colnewyes[colkey]
        print 'LB', lowerbottom[key], 'HT', highertop[key] 
        print 'LT', lowertop[key], 'HB', higherbottom[key]
  
    #return result
    return output

class TestBoxes(unittest.TestCase):
    def setUp(self):
        shape=zeros((9,4,4,4),'f')

        shape[0,0,0,0]=1 #start single cell
        shape[1,:2,0,0]=1 # vent 1 cell top
        shape[2,:2,:2,0]=1 # row entrain 2 cell
        shape[3,:2,:2,:2]=1 # col entrain 4 cell
        shape[4,1:3,:2,:2]=1 # vent/vdet
        shape[5,1:3,1:3,:2]=1 # hent/hdet row
        shape[6,1:3,1:3,1:3]=1 # hent/hdet col
        shape[7,0:4,0:4,0:4]=1 # vent/hent top and bottom  all axes
        shape[8,2,2,2]=1 # vdet/hdet top and bottom  all axes
        self.shape=shape
        self.bxs=boxes(self.shape)
    
    def testTime(self):
        self.assertEqual(self.bxs.shape[0],8)

    def testVentTop(self):
        bxs=self.bxs
        timeslice=bxs[0,...]
        self.assertEqual(timeslice.sum(),2.)
        self.assert_(timeslice[0,0,0,box_id.NOCHG])
        self.assert_(timeslice[1,0,0,box_id.VENT])

    def testHentCol(self):
        bxs=self.bxs
        timeslice=bxs[1,...]
        self.assertEqual(timeslice.sum(),4.)
        self.assert_(timeslice[0,0,0,box_id.NOCHG])
        self.assert_(timeslice[1,0,0,box_id.NOCHG])
        self.assert_(timeslice[0,1,0,box_id.HENT])
        self.assert_(timeslice[1,1,0,box_id.HENT])

    def testHentRow(self):
        bxs=self.bxs
        timeslice=bxs[2,...]
        self.assertEqual(timeslice.sum(),8.)
        self.assert_(timeslice[0,0,0,box_id.NOCHG])
        self.assert_(timeslice[1,0,0,box_id.NOCHG])
        self.assert_(timeslice[0,1,0,box_id.NOCHG])
        self.assert_(timeslice[1,1,0,box_id.NOCHG])
        self.assert_(timeslice[0,0,1,box_id.HENT])
        self.assert_(timeslice[0,1,1,box_id.HENT])
        self.assert_(timeslice[1,0,1,box_id.HENT])
        self.assert_(timeslice[1,1,1,box_id.HENT])

    def testSimVentVdet(self):
        bxs=self.bxs
        timeslice=bxs[3,...]
        self.assertEqual(timeslice.sum(),12.)
        self.assert_(timeslice[0,0,0,box_id.VDET])
        self.assert_(timeslice[0,1,0,box_id.VDET])
        self.assert_(timeslice[0,1,0,box_id.VDET])
        self.assert_(timeslice[0,1,1,box_id.VDET])
        self.assert_(timeslice[1,0,0,box_id.NOCHG])
        self.assert_(timeslice[1,0,1,box_id.NOCHG])
        self.assert_(timeslice[1,1,0,box_id.NOCHG])
        self.assert_(timeslice[1,1,1,box_id.NOCHG])
        self.assert_(timeslice[2,0,0,box_id.VENT])
        self.assert_(timeslice[2,0,1,box_id.VENT])
        self.assert_(timeslice[2,1,0,box_id.VENT])
        self.assert_(timeslice[2,1,1,box_id.VENT])
        
    
    def testSimHentHdetRow(self):
        bxs=self.bxs
        timeslice=bxs[4,...]
        self.assertEqual(timeslice.sum(),12.)
        self.assert_(timeslice[1,0,0,box_id.HDET])
        self.assert_(timeslice[1,0,1,box_id.HDET])
        self.assert_(timeslice[1,1,0,box_id.NOCHG])
        self.assert_(timeslice[1,1,1,box_id.NOCHG])
        self.assert_(timeslice[2,0,0,box_id.HDET])
        self.assert_(timeslice[2,0,1,box_id.HDET])
        self.assert_(timeslice[2,1,0,box_id.NOCHG])
        self.assert_(timeslice[2,1,1,box_id.NOCHG])
        self.assert_(timeslice[1,2,0,box_id.HENT])
        self.assert_(timeslice[1,2,1,box_id.HENT])
        self.assert_(timeslice[2,2,0,box_id.HENT])
        self.assert_(timeslice[2,2,1,box_id.HENT])
        
    def testSimHentHdetCol(self):
        bxs=self.bxs
        timeslice=bxs[5,...]
        self.assertEqual(timeslice.sum(),12.)
        self.assert_(timeslice[1,1,0,box_id.HDET])
        self.assert_(timeslice[1,1,2,box_id.HENT])
        self.assert_(timeslice[1,1,1,box_id.NOCHG])
        self.assert_(timeslice[2,1,0,box_id.HDET])
        self.assert_(timeslice[2,1,2,box_id.HENT])
        self.assert_(timeslice[2,1,1,box_id.NOCHG])
        self.assert_(timeslice[1,2,0,box_id.HDET])
        self.assert_(timeslice[1,2,2,box_id.HENT])
        self.assert_(timeslice[1,2,1,box_id.NOCHG])
        self.assert_(timeslice[2,2,0,box_id.HDET])
        self.assert_(timeslice[2,2,2,box_id.HENT])
        self.assert_(timeslice[2,2,1,box_id.NOCHG])
    
    def testSimEntAllAxes(self):
        bxs=self.bxs
        timeslice=bxs[6,...]
        self.assertEqual(timeslice.sum(),64.)
        self.assertEqual(timeslice[1:3,1:3,1:3,box_id.NOCHG].sum(),8)
        self.assertEqual(timeslice[3,1:3,1:3,box_id.VENT].sum(),4)
        self.assertEqual(timeslice[0,1:3,1:3,box_id.VENT].sum(),4)
        self.assertEqual(timeslice[0:4,0,0:4,box_id.HENT].sum(),16)
        self.assertEqual(timeslice[0:4,3,0:4,box_id.HENT].sum(),16)
        self.assertEqual(timeslice[0:4,1:3,0,box_id.HENT].sum(),8)
        self.assertEqual(timeslice[0:4,1:3,3,box_id.HENT].sum(),8)
        
    def testSimDetAllAxes(self):
        bxs=self.bxs
        timeslice=bxs[7,...]
        self.assertEqual(timeslice.sum(),64.)
        self.assertEqual(timeslice[2,2,2,box_id.NOCHG].sum(),1)
        self.assertEqual(timeslice[0:4,0:4,0:2,box_id.HDET].sum(),32)
        self.assertEqual(timeslice[0:4,0:2,2:4,box_id.HDET].sum(),16)
        self.assertEqual(timeslice[0:4,2:4,3,box_id.HDET].sum(),8)
        self.assertEqual(timeslice[0:4,3,2,box_id.HDET].sum(),4)
        self.assertEqual(timeslice[0:2,2,2,box_id.VDET].sum(),2)
        self.assertEqual(timeslice[3,2,2,box_id.VDET].sum(),1)
        
    def runTest(self):
        pass
        
if __name__=='__main__':
    unittest.main()
