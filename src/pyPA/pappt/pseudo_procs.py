from warnings import warn

from numpy.ma import masked_where, masked_invalid

from pyPA.pappt.lagrangian import boxes, box_id
from PseudoNetCDF import PseudoNetCDFVariable

def simple_pseudo_procs(pa_master, outputfile, spcs, initial, bxs, norm_bxs, contributions, reduce_space):
    #  Assumed order
    #  Pseudo Process: (1) V Detrain (2) V Entrain and Dilution (3) H Detrain (4) H Entrain and Dilution
    # 
    #  Equation subscripts u: nochg; vd: vertical detrainment; ve: vertical entrain and dilution; hd: horizontal detrainment; he: horizontal entrain and dilution
    #  1) (i_u+i_hd+i_vd)/(a_u+a_hd+a_vd)->(i_u+i_hd)/(a_u+a_hd)
    #  2) (i_u+i_hd)/(a_u+a_hd)->(i_u+i_hd+i_ve)/(a_u+a_hd+a_ve)
    #  3) (i_u+i_hd+i_ve)/(a_u+a_hd+a_ve)->(i_u+i_ve)/(a_u+a_ve)
    #  4) (i_u+i_ve)/(a_u+a_ve)->(i_u+i_ve+i_he)/(a_u+a_ve+a_he)
    for spc in spcs:
        init_key = '%s_%s' % (initial, spc)
        if not pa_master.variables.has_key(init_key):
            warn("Missing initial concentration for %s" % spc)
            continue
        var = pa_master.variables[init_key]
        unit = var.units.strip()
        if unit not in contributions:
            warn("Ignoring %s; cannot process unit (%s)" % (spc, unit))
            warn("To add processing for a unit, update the contributions and/or normalizations dictionary with appropriate variable")
            continue

        i=lambda box: reduce_space(masked_where(bxs[..., box] == False, var[:]*contributions[unit])).filled(0)
        a=lambda box: norm_bxs[unit][:,box]
    
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
        hdetrain=masked_invalid(hdetrain).filled(0)
        vdetrain=masked_invalid(vdetrain).filled(0)
        ventrain=masked_invalid(ventrain).filled(0)
        hentrain=masked_invalid(hentrain).filled(0)
    
        del i,a
        
        dilution=vdilution+hdilution
        pseudo_key = 'HENT_%s' % spc
        outputfile.variables[pseudo_key] = PseudoNetCDFVariable(outputfile, pseudo_key, 'f', ('TSTEP',), values = hentrain, units = unit, long_name = 'Horizontally entrained', var_desc = 'Horizontally entrained')
        pseudo_key = 'HDET_%s' % spc
        outputfile.variables[pseudo_key] = PseudoNetCDFVariable(outputfile, pseudo_key, 'f', ('TSTEP',), values = hdetrain, units = unit, long_name = 'Horizontally detrained', var_desc = 'Horizontally detrained')
        pseudo_key = 'VENT_%s' % spc
        outputfile.variables[pseudo_key] = PseudoNetCDFVariable(outputfile, pseudo_key, 'f', ('TSTEP',), values = ventrain, units = unit, long_name = 'Vertically entrained', var_desc = 'Vertically entrained')
        pseudo_key = 'VDET_%s' % spc
        outputfile.variables[pseudo_key] = PseudoNetCDFVariable(outputfile, pseudo_key, 'f', ('TSTEP',), values = vdetrain, units = unit, long_name = 'Vertically detrained', var_desc = 'Vertically detrained')
        pseudo_key = 'EDHDIL_%s' % spc
        outputfile.variables[pseudo_key] = PseudoNetCDFVariable(outputfile, pseudo_key, 'f', ('TSTEP',), values = hdilution, units = unit, long_name = 'Horizontal dilution', var_desc = 'Horizontal dilution')
        pseudo_key = 'EDVDIL_%s' % spc
        outputfile.variables[pseudo_key] = PseudoNetCDFVariable(outputfile, pseudo_key, 'f', ('TSTEP',), values = vdilution, units = unit, long_name = 'Vertical dilution', var_desc = 'Vertical dilution')
    
    
