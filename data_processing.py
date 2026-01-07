from imports import *
from constants import *

def to_radial(pos,centre):
    rel_pos=pos-centre
    radii=np.sqrt(rel_pos[:,0]**2+rel_pos[:,1]**2+rel_pos[:,2]**2)
    return radii

def to_relative(pos,centre):
    rel_pos=pos-centre
    rel_pos_kpc=rel_pos.to(units.kpc)
    return rel_pos_kpc

def get_projected_densities(pos,mass,bin_num,matter_type):

    types=['gas','dm','stars']
    factors=[mass_h_atom,1,1]

    factor=factors[types.index(matter_type)]
    
    xy_binned_mass=stats.binned_statistic_2d(pos[:,0].to_value(units.cm),pos[:,1].to_value(units.cm),mass.to_value(units.g)/factor,bins=[bin_num,bin_num],statistic='sum')
    xz_binned_mass=stats.binned_statistic_2d(pos[:,0].to_value(units.cm),pos[:,2].to_value(units.cm),mass.to_value(units.g)/factor,bins=[bin_num,bin_num],statistic='sum')
    yz_binned_mass=stats.binned_statistic_2d(pos[:,1].to_value(units.cm),pos[:,2].to_value(units.cm),mass.to_value(units.g)/factor,bins=[bin_num,bin_num],statistic='sum')
    
    extent_x=np.max(pos[:,0])-np.min(pos[:,0])
    extent_y=np.max(pos[:,1])-np.min(pos[:,1])
    extent_z=np.max(pos[:,2])-np.min(pos[:,2])

    binwidth_x=extent_x/bin_num
    binwidth_y=extent_y/bin_num
    binwidth_z=extent_z/bin_num

    binarea_xy=binwidth_x.to_value(units.cm)*binwidth_y.to_value(units.cm)
    binarea_xz=binwidth_x.to_value(units.cm)*binwidth_z.to_value(units.cm)
    binarea_yz=binwidth_y.to_value(units.cm)*binwidth_z.to_value(units.cm)

    return xy_binned_mass.statistic/binarea_xy,xz_binned_mass.statistic/binarea_xz,yz_binned_mass.statistic/binarea_yz
