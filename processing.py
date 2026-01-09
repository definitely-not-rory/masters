from imports import *

def to_rel(pos,origin):
    rel=pos-origin
    return rel

def to_rad(pos,origin):
    rel_pos=to_rel(pos,origin)
    radii=np.sqrt(rel_pos[:,0]**2+rel_pos[:,1]**2+rel_pos[:,2]**2)
    return radii

def get_extent(pos):
    dimensions=len(pos[0])
    indexes=np.linspace(0,dimensions-1,dimensions)
    dimension_labels=['x','y','z']
    extents={}
    for index in indexes:
        dimension_min=np.min(pos[:,int(index)])
        dimension_max=np.max(pos[:,int(index)])
        dimension_range=dimension_max-dimension_min
        dimension_extent={'min':dimension_min,'max':dimension_max,'range':dimension_range}
        extents[dimension_labels[int(index)]]=dimension_extent
    return extents
    
