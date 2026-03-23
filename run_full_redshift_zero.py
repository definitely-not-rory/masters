from imports import *
import data_handling as data
from halo_readers import get_snap_num

def run_halo(halo,bin_num):
    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}')
    
    if os.path.isdir(f'figures/{halo}/')!=True:
        os.makedirs(f'figures/{halo}/',exist_ok=True)
    
    figures=['dz_snapshot','fof_scatter',f'{bin_num}/proj_mass_density','radial_mass_density',f'{bin_num}/proj_gas_densities',f'{bin_num}/weighted_mean_gz',f'{bin_num}/nH_col_gz_scatter']
    
    for figure in figures:
        os.makedirs(f'figures/{halo}/{figure}',exist_ok=True)


    data.get_gas_only_data(halo,redshift=0,all_plots=True)

run_halo('T1_Aug',512)
run_halo('T4_Aug',512)
run_halo('halo8',512)