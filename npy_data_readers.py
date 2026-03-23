from imports import *
from halo_readers import get_snap_num

def generate_units_file():
    params=['pos','mass','gz','gmet','nh','rel_pos','radii']
    param_units=[units.Mpc,10**10*units.M_sun,1,1,1,units.Mpc,units.Mpc]
    data=np.array([params,param_units],dtype='object')
    np.save('/cosma/apps/durham/dc-coll7/all_units.npy',data)

def load_units(data,param):
    if os.path.exists('/cosma/apps/durham/dc-coll7/all_units.npy')!=True:
        generate_units_file()
        
    units_data=np.load('/cosma/apps/durham/dc-coll7/all_units.npy',allow_pickle=True)
    all_params=units_data[0]
    all_units=units_data[1]
    if param in all_params:
        param_index=all_params.tolist().index(param)
        param_units=all_units[param_index]
        read_data=data*param_units
    else:
         read_data=data
    return read_data

def read_raw_file(halo,matter_type,param,**kwargs):
    if 'save_dir' in kwargs:
        save_dir=kwargs['save_dir']
    else:
        save_dir='apps'

    if 'snap_num' in kwargs:
        snap_num=str(kwargs['snap_num'])
    elif 'redshift' in kwargs:
        target_redshift=str(kwargs['redshift'])
        snap_num=get_snap_num(halo,target_redshift)[0]
    else:
        sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')

    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
            snap_num='0'+snap_num
    
    raw_data=np.load(f'/cosma/{save_dir}/durham/dc-coll7/halos/{halo}/{snap_num}/raw/{matter_type}/{param}.npy')
    units_data=load_units(raw_data,param)
    return(units_data)

def read_subfind_params(halo,**kwargs):
    if 'save_dir' in kwargs:
        save_dir=kwargs['save_dir']
    else:
        save_dir='apps'
    
    if 'snap_num' in kwargs:
        snap_num=str(kwargs['snap_num'])
    elif 'redshift' in kwargs:
        target_redshift=kwargs['redshift']
        snap_num=str(get_snap_num(halo,target_redshift)[0])
    else:
        sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
            snap_num='0'+snap_num
    
    raw_data=np.load(f'/cosma/{save_dir}/durham/dc-coll7/halos/{halo}/{snap_num}/subfind/halo_params.npy',allow_pickle=True)
    param_names=['redshift','halo_pos','halo_mass','halo_r200','halo_gas_met','halo_star_met','halo_m200']
    param_units=[1,units.Mpc,10**10*units.M_sun,units.Mpc,1,1,10**10*units.M_sun]
    subfind_params={param:raw_data[param_names.index(param)]*param_units[param_names.index(param)] for param in param_names}
    return subfind_params

    