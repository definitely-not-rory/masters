from imports import *

def to_rel(pos,origin):
    rel=pos-origin
    return rel

def to_rad(pos,origin):
    rel_pos=to_rel(pos,origin)
    radii=np.sqrt(rel_pos[:,0]**2+rel_pos[:,1]**2+rel_pos[:,2]**2)
    return radii