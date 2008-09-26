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

from numpy import logical_not, logical_or,logical_and,zeros,where,ones,newaxis,indices,array,arange
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
    done=ones(peak.shape,dtype='int')*(mintop-1)
    
    for lay in range(1,kvvar.shape[kaxis]):
        # for each layer, set the slice to include that layer
        # and all values from times, rows and cols
        kindex[kaxis]=lay-1
        kslice=kvvar[kindex]
        
        # Check if each cell has risen already or if this layer value
        # represents having risen
        critk=pct_of_peak*peak
        risen[...] = logical_or(risen,logical_and(logical_and((kslice>(critk)),lay<4),critk>0))
        
        # If the layer is above the minimum provided layer
        # then check if the value is done
        if lay>=mintop:
            # done is the layer where the value has already 
            # risen, the new value is below the threshold
            # and the layer is greater than the previous done
            done=where(logical_and(logical_or(logical_not(risen),kslice > lower_bound),done==(lay-1)),lay,done)
            
            # This is a straight face check for the algorithm.
            # It is probably no longer necessary
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
    # here is the slice that represents the current layer
    # and below is the slice that represents the layer below
    here=[slice(None)]*4
    below=[slice(None)]*4
    here[kaxis]=slice(1,None)
    below[kaxis]=slice(None,-1)
    
    # dz is the layer thicknesses
    dz=hght.copy()
    dz[here]=dz[here]-hght[below]
    
    # pbl_top will hold the layer index that is the last layer of 
    # thorough mixing and is initialized to the provided minimum
    pbl_top=ones(kvvar.sum(kaxis).shape,dtype='int')*mintop
    
    for lay in range(mintop,kvvar.shape[kaxis]):
        if lay==0:
            pass
        else:
            # here is the slice that represents the current layer
            # below is the slice that represents the layer below
            here[kaxis]=lay
            below[kaxis]=lay-1

            hght_here=dz[here]
            hght_below=dz[below]
            kv=kvvar[below]
            #Crit KV is the vertical diffusivity value of the cell below
            #required to believe that upward transport into this cell
            #signifies that it should be part of the pbl
            #
            # 3% dz_0*avg(dz_0,dz_1)/100
            critkv=0.03*hght_below**2*(1+hght_here/hght_below)/200
            pbl_top=where(logical_and(kv>critkv,pbl_top==(lay-1)),lay,pbl_top)
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
    
    return (where(logical_and(pblhghts[:,:,:,:] <= tops[:,:,:,:], pblhghts[:,:,:,:] > bottoms[:,:,:,:]),1,0)*idx[kaxis]).max(kaxis)

    
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
        self.kv=array( [[12.50309944152832, 48.493377685546875, 86.624908447265625, 58.151683807373047, 2.2295076847076416, 1.2785671949386597, 1.075869083404541, 1.0409003496170044, 1.9714019298553467, 1.1851274967193604, 1.3825576305389404, 1.2397700548171997, 1.0, 1.0, 1.0, 1.0571221113204956], [25.467483520507812, 72.840927124023438, 159.96450805664062, 243.81465148925781, 297.75772094726562, 328.5377197265625, 346.54702758789062, 351.51022338867188, 357.87359619140625, 377.81649780273438, 61.257087707519531, 3.5769078731536865, 2.3378055095672607, 1.0, 1.0, 1.0], [43.479965209960938, 121.53607177734375, 269.95388793945312, 401.28005981445312, 590.94818115234375, 556.8038330078125, 827.30572509765625, 655.28900146484375, 892.768310546875, 710.23773193359375, 751.11358642578125, 668.29766845703125, 10.344052314758301, 4.4217534065246582, 1.5331608057022095, 1.0], [23.319244384765625, 86.194343566894531, 177.90231323242188, 271.3626708984375, 350.075439453125, 387.75286865234375, 450.41873168945312, 437.9488525390625, 485.99362182617188, 470.5758056640625, 474.30386352539062, 475.6044921875, 7.5689301490783691, 1.3417609930038452, 1.0039209127426147, 1.0]] )[:,:,newaxis,newaxis]
        self.z=array( [[33.899448394775391, 84.958000183105469, 170.61279296875, 256.97781372070312, 344.06527709960938, 431.8887939453125, 520.45965576171875, 609.7955322265625, 699.9056396484375, 790.810791015625, 928.682373046875, 1068.42431640625, 1210.0865478515625, 1353.7332763671875, 1597.7078857421875, 1847.6533203125], [33.899448394775391, 84.958000183105469, 170.61279296875, 256.97781372070312, 344.06527709960938, 431.8887939453125, 520.45965576171875, 609.7955322265625, 699.9056396484375, 790.810791015625, 928.682373046875, 1068.42431640625, 1210.0865478515625, 1353.7332763671875, 1597.7078857421875, 1847.6533203125], [33.899448394775391, 84.958000183105469, 170.61279296875, 256.97781372070312, 344.06527709960938, 431.8887939453125, 520.45965576171875, 609.7955322265625, 699.9056396484375, 790.810791015625, 928.682373046875, 1068.42431640625, 1210.0865478515625, 1353.7332763671875, 1597.7078857421875, 1847.6533203125], [33.899448394775391, 84.958000183105469, 170.61279296875, 256.97781372070312, 344.06527709960938, 431.8887939453125, 520.45965576171875, 609.7955322265625, 699.9056396484375, 790.810791015625, 928.682373046875, 1068.42431640625, 1210.0865478515625, 1353.7332763671875, 1597.7078857421875, 1847.6533203125]] )[:,:,newaxis,newaxis]
        self.tops_vertcamx=array( [5, 11, 13, 13] )[:,newaxis,newaxis]
        self.tops_pbldiag=array( [5, 13, 14, 13] )[:,newaxis,newaxis]
        self.hghts_vertcamx=array([431.8887939453125, 1068.42431640625, 1353.7332763671875, 1353.7332763671875])[:,newaxis,newaxis]
        self.hghts_pbldiag=array([431.8887939453125, 1353.7332763671875, 1597.7078857421875, 1353.7332763671875])[:,newaxis,newaxis]


    def testPblDiag(self):
        """
        tops_pbldiag was manually checked to be the right answer
        """
        self.assert_((pbldiag(self.kv,kaxis=1,mintop=5)==self.tops_pbldiag).all())

    def testVertCAMx(self):
        """
        tops_vertcamx was manually checked to be the right answer
        """
        self.assert_((self.tops_vertcamx==vertcamx(self.kv,self.z,kaxis=1,mintop=5)).all())

    def testTops2Shape(self):
        """
        By definition, each horizontal cell should have the top index (0-based)
        pluse 1 layers turned on
        """
        shape=list(self.tops_vertcamx.shape)
        shape[0:1]=[shape[0],self.tops_vertcamx.max()+1]
        self.assert_((tops2shape(self.tops_vertcamx, shape).sum(1)==(self.tops_vertcamx+1)).all())

        shape=list(self.tops_pbldiag.shape)
        shape[0:1]=[shape[0],self.tops_pbldiag.max()+1]
        self.assert_((tops2shape(self.tops_pbldiag, shape).sum(1)==(self.tops_pbldiag+1)).all())


    def testTops2Hghts(self):
        """
        hghts_vertcamx and hghts_pbldiag were manually checked
        """
        self.assert_((tops2hghts(self.tops_vertcamx, self.z,kaxis=1)==self.hghts_vertcamx).all())
        self.assert_((tops2hghts(self.tops_pbldiag, self.z,kaxis=1)==self.hghts_pbldiag).all())
        
    def testPblHghts2Tops(self):
        """
        hghts and tops were both manually checked.
        """
        results=pblhghts2tops(self.hghts_pbldiag[:,newaxis,:,:],self.z,kaxis=1)
        ans=self.tops_pbldiag
        self.assert_((results==ans).all())
        results=pblhghts2tops(self.hghts_vertcamx[:,newaxis,:,:],self.z,kaxis=1)
        ans=self.tops_vertcamx
        self.assert_((results==ans).all())
    
    def runTest(self):
        pass
        
if __name__=='__main__':
    unittest.main()
