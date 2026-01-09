from imports import *
from halo_readers import get_snap_num

def generate_units_file():
    params=['pos','mass','gz','gmet','nh','rel_pos','radii']
    param_units=[units.Mpc,10**10*units.M_sun,z_sol,1,1,units.Mpc,units.Mpc]
    data=np.array([params,param_units],dtype='object')
    np.save('all_units.npy',data)

def load_units(data,param):
    if os.path.exists('all_units.npy')!=True:
        generate_units_file()
        
    units_data=np.load('all_units.npy',allow_pickle=True)
    all_params=units_data[0]
    all_units=units_data[1]

    param_index=all_params.tolist().index(param)
    param_units=all_units[param_index]
    return data*param_units

def read_raw_file(halo,matter_type,param,**kwargs):
    if 'snap_num' in kwargs:
        snap_num=kwargs['snap_num']
    elif 'redshift' in kwargs:
        target_redshift=kwargs['redshift']
        snap_num=get_snap_num(halo,target_redshift)[0]
    else:
        sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    raw_data=np.load(f'halos/{halo}/{snap_num}/raw/{matter_type}/{param}.npy')
    units_data=load_units(raw_data,param)
    return(units_data)

def read_subfind_params(halo,**kwargs):
    if 'snap_num' in kwargs:
        snap_num=kwargs['snap_num']
    elif 'redshift' in kwargs:
        target_redshift=kwargs['redshift']
        snap_num=get_snap_num(halo,target_redshift)[0]
    else:
        sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    raw_data=np.load(f'halos/{halo}/{snap_num}/subfind/halo_params.npy',allow_pickle=True)
    param_names=['redshift','halo_pos','halo_mass','halo_r200']
    param_units=[1,units.Mpc,10**10*units.M_sun,units.Mpc]
    subfind_params={param:raw_data[param_names.index(param)]*param_units[param_names.index(param)] for param in param_names}
    return subfind_params

    