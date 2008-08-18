HeadURL="$HeadURL$"
ChangeDate = "$LastChangedDate$"
RevisionNum= "$LastChangedRevision$"
ChangedBy  = "$LastChangedBy$"
__version__ = RevisionNum

"""
kvextract provides two methods of calculating the planetary boundary 
level: pbldiag and vertcamx.  Each returns a 3D grid (Time, Row, Col).  
It also provides a method of converting the 3D grid to a 4D grid 
(Time, Lay, Row, Col)
"""

from numpy import logical_or,logical_and,zeros,where,ones,newaxis,indices,array,arange
from warnings import warn
import unittest

def pbldiag(kvvar,kaxis=1,mintop=2):
    """
    pbldiag is based on a process developed by Dr. Vizuete and 
    Dr. Kimura at UT and then refined and automated by Barron Henderson

    kvvar - 4D variable of vertical diffusivity values
    kaxis - 0-based layer dimension of kvvar and hght
    mintop - minimum pbl
    
    returns - 3D variable of pbl layer indices
    """
    kindex=[slice(None)]*4
    
    # lower_bound is the minimum vertical diffusivity
    # value.  Under this, is not in the mixing layer.
    lower_bound=2
    
    # pct_of_peak is the threshhold for vertical diffusivity 
    # that must be reached before layer 4 (index 3) (hard coded below)
    pct_of_peak=0.25
    
    # peak is the maximum vertical diffusivity value
    # for each vertical column at each timestep
    peak=kvvar.max(kaxis)
    
    # risen is true if the vertical diffusivity value
    # reaches pct_of_peak*peak before layer 4 (index 3)
    risen=zeros(peak.shape,dtype='bool')
    
    # done is true if the vertical diffusivity value 
    # drops below lower_bound after having risen.
    # it is initialized to the provided minimum value
    done=zeros(peak.shape,dtype='int')
    done[...]=mintop
    
    for lay in range(1,kvvar.shape[kaxis]):
        # for each layer, set the slice to include that layer
        # and all values from times, rows and cols
        kindex[kaxis]=lay
        kslice=kvvar[kindex]
        
        # Check if each cell has risen already or if this layer value
        # represents having risen
        risen[...] = logical_or(risen,logical_and(logical_and((kslice>(pct_of_peak*peak)),lay<3),pct_of_peak*peak>0))
        
        # If the layer is above the minimum provided layer
        # then check if the value is done
        if lay>mintop:
            #done is the layer where the value has already risen, the new value is below the threshold
            #and the layer is greater than the previous done
            done=where(logical_and(logical_and(risen,(kslice < lower_bound)),(done==mintop)),lay,done)
            
            #This is a check for the algorithm.
            if (done==0).any(): raise KeyError
            
    return done

def vertcamx(kvvar,hght,kaxis=1,mintop=2):
    """
    vertcamx comes from an algorithm provided in vertavg 
    on Environ's support tools for CAMx and was implemented
    in Python by Barron Henderson
    
    kvvar - 4D variable of vertical diffusivity values
    hght - 4D variable of layer top heights
    kaxis - 0-based layer dimension of kvvar and hght
    mintop - minimum pbl
    
    returns - 3D variable of pbl layer indices
    """
    # dz_here is the slice that represents the current layer thickness
    # and dz_below is the slice that represents the layer belows thickness
    dz_here=[slice(None)]*4
    dz_below=[slice(None)]*4
    dz_here[kaxis]=slice(1,None)
    dz_below[kaxis]=slice(None,-1)
    
    # dz is the layer thicknesses
    dz=hght
    dz[dz_here]=dz[dz_here]-hght[dz_below]
    
    # kindex_here is the slice that represents the current layer
    # kindex_below is the slice that represents the layer below
    kindex_here=[slice(None)]*4
    kindex_below=[slice(None)]*4
    
    # pbl_top will hold the layer index that is the last layer of 
    # thorough mixing and is initialized to the provided minimum
    pbl_top=zeros(kvvar.sum(kaxis).shape,dtype='int')
    pbl_top[...]=mintop
    
    # notyet is a flag for each vertical column (TIME,ROW,COL)
    # that tells if the column has not yet gone below its minimum
    # vertical diffusivity threshold
    notyet=ones(pbl_top.shape,dtype='bool')
    
    for lay in range(mintop,kvvar.shape[kaxis]):
        if lay==0:
            pass
        else:
            # kindex_here is the slice that represents the current layer
            # kindex_below is the slice that represents the layer below
            kindex_here[kaxis]=lay
            kindex_below[kaxis]=lay-1

            hght_here=dz[kindex_here]
            hght_below=dz[kindex_below]
            kv=kvvar[kindex_below]
            #Crit KV is the vertical diffusivity value of the cell below
            #required to believe that upward transport into this cell
            #signifies that it should be part of the pbl
            #
            critkv=0.03*hght_below**2*(1+hght_here/hght_below)/200
            pbl_top=where(logical_and(kv>critkv,notyet),lay,pbl_top)
            notyet=kv<critkv
    return pbl_top
    
def tops2shape(tops,shape,kaxis=1):
    """
    tops - 3D variable of maximum layers (e.g. TSTEP,ROW,COL)
    shape - iterable of 4 dimensions (e.g. 24,28,65,83)
    kaxis - 0-based layer dimension of shape (e.g. 1)
    """
    # idx(TLRC,TSTEP,LAY,ROW,COL)
    idx=indices(shape)
    
    # Create a slice for indexing
    tops_idx=[slice(None),slice(None),slice(None),slice(None)]
    tops_idx[kaxis]=newaxis
    
    # Where tops is less than (or equal) to the layer in idx
    # Return true
    return idx[kaxis]<=tops[tops_idx]

def pblhghts2tops(pblhghts,hghts,kaxis=1):
    """
    pblhghts - 4D variable of maximum layers (e.g. TSTEP,LAY=1,ROW,COL)
    hght - 4D variable of layer top heights
    kaxis - 0-based layer dimension of hghts
    """
    tops=hghts
    bottoms=zeros(hghts.shape,'f')
    bottoms[:,1:,:,:]=hghts[:,:-1,:,:]

    idx=indices(hghts.shape)
    
    return (where(logical_and(pblhghts<=tops,pblhghts>bottoms),1,0)*idx[kaxis]).max(kaxis)

    
def tops2hghts(tops,hghts,kaxis=1):
    """
    tops - 3D variable of maximum layers (e.g. TSTEP,ROW,COL)
    hght - 4D variable of layer top heights
    kaxis - 0-based layer dimension of hghts
    """
    # 1) Convert tops to shape and multiply by heights
    # 2) take max height from each Time, Row, Col
    return (tops2shape(tops,hghts.shape,kaxis)*hghts).max(kaxis)

class TestPbl(unittest.TestCase):
    def setUp(self):
        warn("Test case is not fully implemented")
        self.pblhghts=array(arange(1,5,dtype='f').reshape((2,2))*10.,ndmin=3).repeat(4,0)+5
        self.hghts=array(arange(1,5,dtype='f')*20.,ndmin=4).swapaxes(1,3).repeat(2,3).repeat(2,2).repeat(4,0)

    def testPblDiag(self):
        pass

    def testVertCAMx(self):
        pass

    def testTops2Shape(self):
        pass

    def testTops2Hghts(self):
        pass
        
    def testPblHghts2Tops(self):
        results=pblhghts2tops(self.pblhghts,self.hghts)
        ans=array([[[0, 1],[1, 2]]]*4)
        self.assert_((results==ans).all())
    
    def runTest(self):
        pass
if __name__=='__main__':
    unittest.main()