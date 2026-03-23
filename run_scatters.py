from imports import *
import plot_generation as plot
from halo_readers import get_snap_num

def run_scatters(halo,bin_num,start_redshift,end_redshift):
    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}')
    
    start_snap_num=int(get_snap_num(halo,start_redshift)[0])
    end_snap_num=int(get_snap_num(halo,end_redshift)[0])
    
    snapshots=np.arange(start_snap_num,end_snap_num+1)

    for snapshot in snapshots:
        plot.nH_col_gz_scatter(halo,512,'xy',snap_num=snapshot)
        plot.rho_gz_scatter(halo,512,'xy',snap_num=snapshot)

