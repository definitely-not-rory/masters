from imports import *
import data_handling as data
from halo_readers import get_snap_num

def run_halo(halo,bin_num,start_redshift,end_redshift):
    start_snap_num=int(get_snap_num(halo,start_redshift)[0])
    end_snap_num=int(get_snap_num(halo,end_redshift)[0])
    
    snapshots=np.arange(start_snap_num,end_snap_num+1)

    if os.path.isdir(f'figures/{halo}/{bin_num}px/')!=True:
        os.makedirs(f'figures/{halo}/{bin_num}px/',exist_ok=True)