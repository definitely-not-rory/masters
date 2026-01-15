from imports import *
import data_handling as data
from halo_readers import get_snap_num

def run_halo(halo,bin_num,start_redshift,end_redshift):
    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}')
    
    start_snap_num=int(get_snap_num(halo,start_redshift)[0])
    end_snap_num=int(get_snap_num(halo,end_redshift)[0])
    
    snapshots=np.arange(start_snap_num,end_snap_num+1)

    if os.path.isdir(f'figures/{halo}/')!=True:
        os.makedirs(f'figures/{halo}/',exist_ok=True)
    
    figures=['dz_snapshot','fof_scatter',f'{bin_num}/proj_mass_density','radial_mass_density',f'{bin_num}/proj_gas_densities',f'{bin_num}/weighted_mean_gz',f'{bin_num}/nH_col_gz_scatter']
    
    for figure in figures:
        os.makedirs(f'figures/{halo}/{figure}',exist_ok=True)

    for snapshot in snapshots:
        data.get_gas_only_data(halo,snap_num=snapshot)
    
    data.get_threshold_behaviour_data(halo,bin_num=bin_num,redshifts=[start_redshift,end_redshift])

run_halo('halo8',512,4,1)