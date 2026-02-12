from imports import *
from halo_readers import get_halos, get_redshifts, get_snap_num, get_redshift, get_snap_nums
from npy_data_readers import read_raw_file, read_subfind_params
import plot_generation as plot 
import processing as calc  

def get_raw_data(halo,**kwargs): #Function to import all raw subfind and snapshot data into .npy files for a given halo and snapshot number/redshift
    
    #--- Halo Selection ---

    if 'path' in kwargs: #Detects for provision of alternate data directory
        alt_path=kwargs['path']
        available_halos=get_halos(path=alt_path) #Loads halos in alternate directory
    else:
        available_halos=get_halos() #Loads halos in default directory if no alternate is provided
    
    if halo not in available_halos: #Detects if requested halo exists
        print('Halo does not exist, please select from:') #Returns lis of available halo in provided/default directory
        print(*available_halos,sep='\n')
        sys.exit('Halo Name Not Found in Directory') #Error message for incorrect halo name
    
    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}') #Create halo directory if required
        print(f'Created directory for {halo}')
    else:
        print(f'\nDirectory for {halo} located')

    #--- Snapshot Number Selection ---
    get_redshifts(halo)
    if 'snap_num' in kwargs: #Detects if a specific snapshot number is provided
        snap_num=str(kwargs['snap_num']) #Loads snapshot number integer input
        while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
            snap_num='0'+snap_num
        print(f'\nSnapshot Number {snap_num} Provided')
    elif 'redshift' in kwargs: #Detects if a target redshift is provided
        target_redshift=np.float64(kwargs['redshift']) #Loads provided target redshift
        snap_nums, snap_redshifts=get_redshifts(halo)
        snap_redshifts=np.array(snap_redshifts).astype(np.float64)
        if target_redshift==0: #If requested redshift is z=0, use last possible snapshot
            snap_num=snap_nums[-1]
            snap_redshift=np.float64(snap_redshifts[-1])
            print(f'\nTarget Redshift: 0, Using Final Snapshot ({snap_num})')
        else:
            snap_num,snap_redshift=get_snap_num(halo,target_redshift)
    else:
        sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")') #Error message for incorrect snapshot number/redshift

    if 'plot_dz_snap' in kwargs: #Detects if plots are enabled
        if kwargs['plot_dz_snap']==True:
            if os.path.isdir(f'figures/{halo}/dz_snapshot')!=True:
                os.makedirs(f'figures/{halo}/dz_snapshot',exist_ok=True)
            if 'redshift' in kwargs:
                plot.dz_snapshot(halo,snap_num,target_redshift=target_redshift,snap_redshift=snap_redshift) #Generates dz vs snapshot plot with target redshift elements if redshift was used instead of snapshot number
            else:
                plot.dz_snapshot(halo,snap_num) #Generates dz vs snapshot plot without target redshift elements if snapshot number was used
    elif 'all_plots' in kwargs:
        if os.path.isdir(f'figures/{halo}/dz_snapshot')!=True:
                os.makedirs(f'figures/{halo}/dz_snapshot',exist_ok=True)
        if kwargs['all_plots']==True: #Checks if override for all plots is enabled
            if 'redshift' in kwargs:
                plot.dz_snapshot(halo,snap_num,target_redshift=target_redshift,snap_redshift=snap_redshift) #Generates dz vs snapshot plot with target redshift elements if redshift was used instead of snapshot number
            else:
                plot.dz_snapshot(halo,snap_num) #Generates dz vs snapshot plot without target redshift elements if snapshot number was used


    if 'path' in kwargs: #Detects for provision of alternate data directory
        loc=kwargs['path']
    else:
        loc='/cosma8/data/dp004/lyra/original_sample/' #Uses default location if none provided
    suffix='/output/' #Suffix to ensure data from each snapshot is read

    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}')!=True: #Checks for dedicated directory for chosen snapshot
        print(f'\nNo Directory Detected for {halo} Snapshot {snap_num}') 
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}') #Creates snapshot directory if necessary
        print(f'{halo} Snapshot {snap_num} Directory Created')
    else:
        print(f'\n{halo} Snapshot {snap_num} Directory Located')    

    #--- Subfind Data Handling ---

    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/subfind')!=True: #Checks for subfind data directory within snapshot directory
        print(f'\nNo {halo} Snapshot {snap_num} subfind directory detected')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/subfind') #Creates subfind directory if necessary
        print(f'Subfind directory for {halo} Snapshot {snap_num} created')
    else:
        print(f'\n{halo} Snapshot {snap_num} Subfind Directory Located')

    print(f'\nImporting {halo} Snapshot {snap_num} Subfind Data')
    subfind_data = ar.gadget_subfind.load_subfind(int(snap_num), dir=loc + halo + suffix) #Import subfind dataset
    print('Subfind Data imported')

    if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/subfind/fof_positions.npy')!=True or os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/subfind/fof_masses.npy')!=True or os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/subfind/halo_params.npy')!=True: #Detects if subfind data arrays have been already imported and stored external
        print(f'\nNo Subfind FoF Data Located\Loading subfind data for {halo} from Snapshot {snap_num}') #Creates relevant external data arrays if required

        #--- Subfind-Wide Data Processing ---
        sf_positions=subfind_data.data['fpos']  #FoF Positions (in Mpc)
        all_sf_masses=subfind_data.data['fmty'] #FoF Masses Across All Types (in 10^10 M_sol)

        high_res_mask=(all_sf_masses[:,2]==0) & (all_sf_masses[:,3]==0) #Masks mass data to exclude low-res structures
        high_res_sf_masses=all_sf_masses[high_res_mask]
        sf_masses=np.sum(high_res_sf_masses,axis=1) #Calculates total high-res mass for each FoF group
        
        np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/subfind/fof_positions.npy',sf_positions) #Saves subfind-wide arrays to external .npy files
        print('Subfind FoF position data saved')

        np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/subfind/fof_masses.npy',sf_masses)
        print('Subfind FoF mass data saved')

        #--- FoF Group 1/Main Halo Data Processing ---
        halo_pos=sf_positions[0] #Main Halo Position (in Mpc)
        halo_mass=sf_masses[0] #Main Halo Mass (in 10^10 M_sol)
        halo_r200=subfind_data.data['frc2'][0] #Main Halo R_200_crit (in Mpc)

        redshift=subfind_data.redshift #Retrieve snapshot redshift value from header
        
        saved_params=np.array([redshift,halo_pos,halo_mass,halo_r200],dtype='object') #Create storage array for FoF Group 1 parameters
        
        np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/subfind/halo_params.npy',saved_params) #Save FoF Group 1 parameters to external .npy file
        print('Subfind Halo Parameters Saved')
    else:
        print('\nFoF Subfind Data Located')
    
    if 'plot_fof_scatter' in kwargs: #Detects if plots are enabled
        if os.path.isdir(f'figures/{halo}/fof_scatter')!=True:
                os.makedirs(f'figures/{halo}/fof_scatter',exist_ok=True)
        if kwargs['plot_fof_scatter']==True:
            plot.fof_scatter(halo,snap_num=snap_num) #Generates FoF scatter plot
    elif 'all_plots' in kwargs: #Checks if override for all plots is enabled
        if os.path.isdir(f'figures/{halo}/fof_scatter')!=True:
                os.makedirs(f'figures/{halo}/fof_scatter',exist_ok=True)
        if kwargs['all_plots']==True:
            plot.fof_scatter(halo,snap_num=snap_num)

    #--- Snapshot Data Handling ---

    #--- Directory Management ---
    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw')!=True: #Checks for raw snapshot data directory within snapshot directory
        print(f'\nNo {halo} raw Snapshot {snap_num} data directory detected')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw') #Creates raw data directory if necessary

        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/gas') #Creates raw matter type directories
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/dm')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/stars')

        print(f'Directory for {halo} raw Snapshot {snap_num} data created')
    else:
        print(f'\n{halo} raw Snapshot {snap_num} data Directory Located')

    #--- Data Import Selection ---
    #Arrays of necessary snapshot data needed to be saved externally
    gas_data=['pos','mass','gz','gmet','nh'] #Gas: Position (Mpc), Mass (10^10 M_sol), Metallicity ( ), Gas Mass Fractions ( ), Neutral Hydrogen Mass Fraction ( )
    dm_data=['pos','dm_params','mass'] #DM: Position (Mpc), Hubble Parameter and DM Particle Mass (10^10 M_sol), Mass (10^10 M_sol)
    stars_data=['pos','mass'] #Stars: Position (Mpc), Mass (10^10 M_sol)
    req_snapshot_data=[gas_data,dm_data,stars_data] #2D array for iterative efficiency
    
    types=['Gas','DM','Stars'] #Text output labels
    type_directories=['gas','dm','stars'] #names of target directories to save .npy files to
    lyra_types=[0,1,4] #Type numbers for gas/DM/stars in raw Lyra data

    current_gas_data=[file[:-4] for file in os.listdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/gas')] #Import list of current files saved for each file
    current_dm_data=[file[:-4] for file in os.listdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/dm')]
    current_stars_data=[file[:-4] for file in os.listdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/stars')]

    current_data=[current_gas_data,current_dm_data,current_stars_data] #2D array for comparison to required data

    #--- Data Imports ---
    for matter_type in req_snapshot_data: #Check whether all required data is present in each type directory
        type_index=req_snapshot_data.index(matter_type) #Determine which index (0,1,2) current target matter type is registered in
        type_name=types[type_index] 
        if set(matter_type).issubset(current_data[type_index])!=True: #Determine if current directory entries matches required data
            lyra_type=lyra_types[type_index] #Import Lyra type for current target type

            #--- Snapshot Importing ---
            print(f'\n{halo} Snapshot {snap_num} {type_name} raw data incomplete\nImporting Snapshot {snap_num} {type_name} data') #If current data is incomplete, import relevant type's snapshot data
            snapshot_data= ar.gadget_snap.gadget_snapshot(loc + halo + suffix + 'snapdir_'+str(snap_num)+'/snapshot_'+str(snap_num),loadonlytype=[lyra_type],loadonlyhalo = 0, subfind = subfind_data, hdf5=True, lazy_load=True)
            print(f'{halo} Snapshot {snap_num} {type_name} data loaded')

            print(f'\nSaving required {type_name} snapshot data to external file')

            #--- Initial Raw Data Handling and External .npy Storage ---
            if type_index==1: #Applies necessary exceptions for Lyra dark matter data
                dm_pos=snapshot_data.data['pos'] #Loads DM position data
                print(f'\nDM pos data loaded')

                np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/dm/pos.npy',dm_pos) #Saves DM position data to external .npy file
                print(f'DM pos data saved externally')

                #Loads necessary parameters to determine constant DM particle mass
                h=snapshot_data.header['HubbleParam'] #Import h and header DM particle mass
                raw_dm_particle_mass=snapshot_data.header['MassTable'][1]
                if raw_dm_particle_mass!=0:
                    dm_particle_mass=raw_dm_particle_mass/h
                else:
                    nonzero_dm_particle_mass=read_raw_file('T1_Aug','dm','dm_params',snap_num=152)[1]
                    dm_particle_mass=nonzero_dm_particle_mass/h #Apply relevant conversion from header DM particle mass to 'real' mass
                print(f'\nDM Mass Parameters data loaded')

                dm_params=np.array([h,dm_particle_mass]) #Stores all DM mass parameters in array for external file
                np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/dm/dm_params.npy',dm_params) #Saves DM parameters to external .npy file
                print(f'DM Mass Parameters data saved externally')

                dm_mass=np.full_like(dm_pos[:,0],dm_particle_mass) #Creates 'normal' DM mass array using params and loaded snapshot data
                print(f'\nDM mass data calculated')

                np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/dm/mass.npy',dm_mass) #Saves DM mass array
                print(f'DM mass data saved externally')     
            else: #Completes all non-exception exports
                to_save=req_snapshot_data[type_index] #Loads all requires snapshot data labels
                for param in to_save: #Exports each required parameter for target type
                    data=snapshot_data.data[param] #Loads target type's target parameter from snapshot
                    print(f'\n{type_name} {param} data loaded')

                    directory=type_directories[type_index] #Finds correct type directory for data storage
                    np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/{directory}/{param}.npy',data) #Saves parameter data to external .npy file
                    print(f'{type_name} {param} data saved externally')
        else:
            print(f'\n{halo} Snapshot {snap_num} {type_name} raw data complete')
    print(f'\n{halo} Snapshot {snap_num} Raw Data Import Complete')

def get_mass_density_data(halo,**kwargs):
    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}') #Create halo directory if required
        print(f'Created directory for {halo}')
    else:
        print(f'\nDirectory for {halo} located')

    if 'snap_num' not in kwargs:
        if 'redshift' in kwargs:
            target_redshift=kwargs['redshift']
            snap_num,snap_redshift=get_snap_num(halo,target_redshift)
        else:
            sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    else:
        snap_num=kwargs['snap_num']
        snap_redshift=get_redshift(halo,snap_num)


    snap_num=str(snap_num)
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
            snap_num='0'+snap_num

    type_directories=['gas','dm','stars']
    params=[['pos','mass'],['pos','mass'],['pos','mass']]
    type_names=['Gas','DM','Stars']
    
    loaded_data = {}

    for matter_type in type_directories:
        type_index=type_directories.index(matter_type)
        type_name=type_names[type_index]
        to_load=params[type_index]

        if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/{matter_type}/pos.npy')!=True or os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/{matter_type}/mass.npy')!=True:
            print(f'\nSnapshot {snap_num} {type_name} Raw Data Incomplete')
            if 'all_plots' in kwargs: #Checks if override for all plots is enabled
                if kwargs['all_plots']==True:
                    get_raw_data(halo,redshift=target_redshift,all_plots=True)
                else:
                    get_raw_data(halo,redshift=target_redshift)
        else:
            print(f'\nSnapshot {snap_num} {type_name} raw data located')
        
        halo_pos=read_subfind_params(halo,snap_num=snap_num)['halo_pos']
        redshift=read_subfind_params(halo,snap_num=snap_num)['redshift']

        type_data={param:read_raw_file(halo,matter_type,param,snap_num=snap_num) for param in to_load}
        
        loaded_data[matter_type]=type_data

    for matter_type, data in loaded_data.items():
        type_name=type_names[type_directories.index(matter_type)]

        if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/{matter_type}/rel_pos.npy')!=True:
            print(f'{type_name} Relative Position Data Not Found, Generating')

            rel_pos=calc.to_rel(data['pos'],halo_pos)
            loaded_data[matter_type]['rel_pos']=rel_pos

            np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/{matter_type}/rel_pos.npy',rel_pos.value)
            print(f'{type_name} Relative Position File Generated')
        else:
            rel_pos=read_raw_file(halo,matter_type,'rel_pos',snap_num=snap_num)
            print(f'{type_name} Relative Positions Loaded')

            loaded_data[matter_type]['rel_pos']=rel_pos

    for matter_type, data in loaded_data.items():
        type_name=type_names[type_directories.index(matter_type)]
        if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/{matter_type}/radii.npy')!=True:
            print(f'{type_name} Radial Position Data Not Found, Generating')

            radii=calc.to_rad(data['pos'],halo_pos)
            loaded_data[matter_type]['radii']=radii

            np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/{matter_type}/radii.npy',radii.value)
            print(f'{type_name} Radial Position File Generated')
        else:
            radii=read_raw_file(halo,matter_type,'radii',snap_num=snap_num)
            print(f'{type_name} Radial Positions Loaded')

            loaded_data[matter_type]['radii']=radii
    
    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned')!=True:
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned')
                      
    if 'bins' in kwargs:
        bin_num=kwargs['bins']
    else:
        bin_num=512

    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px')!=True:
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px')
    
    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/total_mass')!=True:
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/total_mass')
    
    for matter_type, data in loaded_data.items():
        type_name=type_names[type_directories.index(matter_type)]
        if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/total_mass/{matter_type}.npy')!=True:
            extents=calc.get_extent(data['rel_pos'])
            print(f'\n{matter_type} spatial extent calculated')
            
            bin_widths=[extents[dimension]['range'].to_value(units.cm)/bin_num for dimension in extents]
            print(f'\n{matter_type} bin widths calculated')

            plane_indexes=[[0,1],[0,2],[1,2]]
            plane_names=['xy','xz','yz']

            cm_rel_pos=data['rel_pos'].to_value(units.cm)
            g_mass=data['mass'].to_value(units.g)

            proj_densities={}

            for plane in plane_indexes:
                binned_mass=stats.binned_statistic_2d(cm_rel_pos[:,plane[0]],cm_rel_pos[:,plane[1]],g_mass,bins=[bin_num,bin_num],statistic='sum').statistic
                print(f'\n{matter_type} {plane_names[plane_indexes.index(plane)]} binned masses calculated')
                bin_area=bin_widths[plane[0]]*bin_widths[plane[1]]
                
                proj_mass_dens=binned_mass/bin_area
                print(f'\n{matter_type} {plane_names[plane_indexes.index(plane)]} projected densities calculated')

                proj_densities[plane_names[plane_indexes.index(plane)]]=proj_mass_dens

            loaded_data[matter_type]['proj_dens']=proj_densities

            save_data=[proj_densities['xy'],proj_densities['xz'],proj_densities['yz']]
            np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/total_mass/{matter_type}.npy',save_data)
            print(f'\n{matter_type} projected density file created')
        else:
            loaded_data[matter_type]['proj_dens']=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/total_mass/{matter_type}.npy')
            print(f'\n{matter_type} projected density file loaded')

    if 'plot_proj_dens' in kwargs: #Detects if plots are enabled
        if kwargs['plot_proj_dens']==True:
            if os.path.isdir(f'figures/{halo}/{bin_num}/proj_mass_density')!=True:
                os.makedirs(f'figures/{halo}/{bin_num}/proj_mass_density',exist_ok=True)
            plot.proj_mass_density(halo,bin_num,snap_num=snap_num) #Generates g/cm^2 projected density plot
    elif 'all_plots' in kwargs: #Checks if override for all plots is enabled
        if kwargs['all_plots']==True:
            if os.path.isdir(f'figures/{halo}/{bin_num}/proj_mass_density')!=True:
                os.makedirs(f'figures/{halo}/{bin_num}/proj_mass_density',exist_ok=True)
            plot.proj_mass_density(halo,bin_num,snap_num=snap_num)

    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/radial_mass_density')!=True:
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/radial_mass_density')

    for matter_type, data in loaded_data.items():
        type_name=type_names[type_directories.index(matter_type)]
        if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/radial_mass_density/{matter_type}.npy')!=True:
            spherical_bins = np.logspace(-4,-1,50)*units.Mpc
            bin_centres=(spherical_bins[:-1] + spherical_bins[1:])/2

            sphere_volumes=[4/3*np.pi*radii.value**3 for radii in spherical_bins]
            shell_volumes=np.diff(sphere_volumes)
            
            radial_densities=stats.binned_statistic(data['radii'],data['mass'],bins=spherical_bins,statistic='sum').statistic/shell_volumes*units.M_sun/units.Mpc**3

            kpc_bin_centres=bin_centres.to_value(units.kpc)
            gcm3_radial_densities=radial_densities.to_value(units.g/units.cm**3)

            radial_density_profile=np.array([kpc_bin_centres,gcm3_radial_densities])
            np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/radial_mass_density/{matter_type}.npy',radial_density_profile)
    
    if 'plot_rad_dens' in kwargs: #Detects if plots are enabled
        if kwargs['plot_rad_dens']==True:
            if os.path.isdir(f'figures/{halo}/radial_mass_density')!=True:
                os.makedirs(f'figures/{halo}/radial_mass_density',exist_ok=True)
            plot.radial_mass_density(halo,snap_num=snap_num) #Generates g/cm^2 radial density plot
    elif 'all_plots' in kwargs: #Checks if override for all plots is enabled
        if kwargs['all_plots']==True:
            if os.path.isdir(f'figures/{halo}/radial_mass_density')!=True:
                os.makedirs(f'figures/{halo}/radial_mass_density',exist_ok=True)
            plot.radial_mass_density(halo,snap_num=snap_num)

            
def get_gas_only_data(halo,**kwargs):    
    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}') #Create halo directory if required
        print(f'Created directory for {halo}')
    else:
        print(f'\nDirectory for {halo} located')

    if 'snap_num' not in kwargs:
        if 'redshift' in kwargs:
            target_redshift=kwargs['redshift']
            snap_num,snap_redshift=get_snap_num(halo,target_redshift)
        else:
            sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    else:
        snap_num=kwargs['snap_num']
        snap_redshift=get_redshift(halo,snap_num)  

    snap_num=str(snap_num)
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
            snap_num='0'+snap_num


    if 'bins' in kwargs:
        bin_num=kwargs['bins']
    else:
        bin_num=512 

    to_load=['gmet','mass','nh','gz']
    loaded_data={}
    ran_get_raw=False
    
    for param in to_load:
        if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/gas/{param}.npy')!=True:
            print(f'\nSnapshot {snap_num} {param} Raw Data Incomplete')
            if 'all_plots' in kwargs: #Checks if override for all plots is enabled
                if kwargs['all_plots']==True:
                    print('Getting raw data with plots')
                    get_raw_data(halo,bin_num=bin_num,snap_num=snap_num,all_plots=True)
            else:
                print('Getting raw data without plots')
                get_raw_data(halo,bin_num=bin_num,snap_num=snap_num)
        else:
            print(f'\nSnapshot {snap_num} {param} raw data located')
            
        loaded_data[param]=read_raw_file(halo,'gas',param,snap_num=snap_num)
        
    loaded_data['halo_r200']=read_subfind_params(halo,snap_num=snap_num)['halo_r200']  

    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned')!=True:
        if 'all_plots' in kwargs: #Checks if override for all plots is enabled
                if kwargs['all_plots']==True:
                    get_mass_density_data(halo,bin_num=bin_num,snap_num=snap_num,all_plots=True)
        else:
            get_mass_density_data(halo,bin_num=bin_num,snap_num=snap_num)

    if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/raw/gas/rel_pos.npy')!=True:
        get_mass_density_data(halo,bin_num=bin_num,snap_num=snap_num)

    loaded_data['rel_pos']=read_raw_file(halo,'gas','rel_pos',snap_num=snap_num)
    loaded_data['total_mass']=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/total_mass/gas.npy')

    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only')!=True:
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only')

    req_densities=['total_mass','hydrogen_mass','nH_mass','nH_col']
    current_densities=[file[:-4] for file in os.listdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only')]

    if set(req_densities).issubset(current_densities)!=True:
        np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/total_mass.npy',loaded_data['total_mass'])

        extents=calc.get_extent(loaded_data['rel_pos'])
        print(f'\nSpatial extent calculated')
            
        bin_widths=[extents[dimension]['range'].to_value(units.cm)/bin_num for dimension in extents]
        print(f'\nBin widths calculated')

        plane_indexes=[[0,1],[0,2],[1,2]]
        plane_names=['xy','xz','yz']

        cm_rel_pos=loaded_data['rel_pos'].to_value(units.cm)
        g_mass=loaded_data['mass'].to_value(units.g)
        H_mass_frac=loaded_data['gmet'][:,0]
        nH_mass_frac=loaded_data['nh']
        
        proj_densities={'hydrogen_mass':{},'nH_mass':{},'nH_col':{}}

        for plane in plane_indexes:
            plane_name=plane_names[plane_indexes.index(plane)]
            
            binned_hydrogen_mass=stats.binned_statistic_2d(cm_rel_pos[:,plane[0]],cm_rel_pos[:,plane[1]],g_mass*H_mass_frac,bins=[bin_num,bin_num],statistic='sum').statistic
            binned_nH_mass=stats.binned_statistic_2d(cm_rel_pos[:,plane[0]],cm_rel_pos[:,plane[1]],g_mass*H_mass_frac*nH_mass_frac,bins=[bin_num,bin_num],statistic='sum').statistic
            binned_nH_col=stats.binned_statistic_2d(cm_rel_pos[:,plane[0]],cm_rel_pos[:,plane[1]],g_mass*H_mass_frac*nH_mass_frac/mass_H_atom,bins=[bin_num,bin_num],statistic='sum').statistic

            bin_area=bin_widths[plane[0]]*bin_widths[plane[1]]

            hydrogen_mass_density=binned_hydrogen_mass/bin_area
            nH_mass_density=binned_nH_mass/bin_area
            nH_col_density=binned_nH_col/bin_area

            proj_densities['hydrogen_mass'][plane_name]=hydrogen_mass_density
            proj_densities['nH_mass'][plane_name]=nH_mass_density
            proj_densities['nH_col'][plane_name]=nH_col_density
        
        for density in proj_densities:
            save_data=[proj_densities[density]['xy'],proj_densities[density]['xz'],proj_densities[density]['yz']]
            np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/{density}.npy',save_data)
    else:
        print('Gas Density Files Located')

    if 'plot_gas_dens' in kwargs: #Detects if plots are enabled
        if kwargs['plot_gas_dens']==True:
            if os.path.isdir(f'figures/{halo}/{bin_num}/proj_gas_densities')!=True:
                os.makedirs(f'figures/{halo}/{bin_num}/proj_gas_densities',exist_ok=True)
            plot.proj_gas_densities(halo,bin_num,snap_num=snap_num) #Generates gas density plot
    elif 'all_plots' in kwargs: #Checks if override for all plots is enabled
        if kwargs['all_plots']==True:
            if os.path.isdir(f'figures/{halo}/{bin_num}/proj_gas_densities')!=True:
                os.makedirs(f'figures/{halo}/{bin_num}/proj_gas_densities',exist_ok=True)
            plot.proj_gas_densities(halo,bin_num,snap_num=snap_num)

    if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/mean_gz.npy')!=True or os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/bin_radii.npy')!=True or os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/bin_masses.npy')!=True:
        plane_indexes=[[0,1],[0,2],[1,2]]
        plane_names=['xy','xz','yz']

        rel_pos=loaded_data['rel_pos']
        metallicity=loaded_data['gz']/z_sol
        mass=loaded_data['mass']

        total_mass=np.sum(mass)

        mean_metallicties={}

        for plane in plane_indexes:
            plane_name=plane_names[plane_indexes.index(plane)]
            mean_metallicties[plane_name]={}

            binned_masses=stats.binned_statistic_2d(rel_pos[:,plane[0]],rel_pos[:,plane[1]],mass,bins=[bin_num,bin_num],statistic='sum')

            bin_grid=np.meshgrid((binned_masses.x_edge[:-1]+binned_masses.x_edge[1:])/2,(binned_masses.y_edge[:-1]+binned_masses.y_edge[1:])/2)
            bin_radii=np.sqrt(bin_grid[0]**2+bin_grid[1]**2)
            mean_metallicties[plane_name]['bin_radii']=bin_radii

            bin_masses=binned_masses.statistic
            mean_metallicties[plane_name]['bin_masses']=bin_masses

            weighted_metallicity=stats.binned_statistic_2d(rel_pos[:,plane[0]],rel_pos[:,plane[1]],metallicity*mass,bins=[bin_num,bin_num],statistic='sum').statistic
            weighted_mean_metallicity=np.divide(weighted_metallicity,bin_masses)
            mean_metallicties[plane_name]['mean_gz']=weighted_mean_metallicity

        to_save=['mean_gz','bin_radii','bin_masses']
        save_data=np.array([[mean_metallicties['xy'][param],mean_metallicties['xz'][param],mean_metallicties['yz'][param]] for param in to_save])

        np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/mean_gz.npy',save_data[0])
        np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/bin_radii.npy',save_data[1])
        np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/bin_masses.npy',save_data[2])

    if 'plot_mean_gz' in kwargs: #Detects if plots are enabled
        if kwargs['plot_mean_gz']==True:
            if os.path.isdir(f'figures/{halo}/{bin_num}/weighted_mean_gz')!=True:
                os.makedirs(f'figures/{halo}/{bin_num}/weighted_mean_gz',exist_ok=True)
            plot.weighted_mean_gz(halo,bin_num,snap_num=snap_num) #Generates gas density plot
    elif 'all_plots' in kwargs: #Checks if override for all plots is enabled
        if kwargs['all_plots']==True:
            if os.path.isdir(f'figures/{halo}/{bin_num}/weighted_mean_gz')!=True:
                os.makedirs(f'figures/{halo}/{bin_num}/weighted_mean_gz',exist_ok=True)
            plot.weighted_mean_gz(halo,bin_num,snap_num=snap_num)

    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked')!=True:
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked')   
    
    obj_types=['DLA','subDLA','LymanLimit','lo_z_DLA','lo_z_subDLA','lo_z_LymanLimit','lo_z']
    current_masks=[file[:-4] for file in os.listdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked')]

    if set(obj_types).issubset(current_masks)!=True:
        nH_col=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/nH_col.npy')
        mean_gz=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/mean_gz.npy')
        
        masks={'DLA':np.log10(nH_col) >20.3,'subDLA': (20.3>np.log10(nH_col)) & (np.log10(nH_col) >19),'LymanLimit': (19>np.log10(nH_col)) & (np.log10(nH_col) >17.2)}
        
        lo_z_mask=np.log10(mean_gz)<-3
        lo_z_nH_col=np.ma.masked_where(~lo_z_mask,nH_col)
        lo_z_mean_gz=np.ma.masked_where(~lo_z_mask,mean_gz)
        
        np.savez(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked/lo_z.npz', nH_col_data=lo_z_nH_col.data, nH_col_mask=lo_z_nH_col.mask, mean_gz_data=lo_z_mean_gz.data, mean_gz_mask=lo_z_mean_gz.mask)


        for mask in masks:
            masked_nH_col=np.ma.masked_where(~masks[mask],nH_col)
            masked_mean_gz=np.ma.masked_where(~masks[mask],mean_gz)

            np.savez(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked/{mask}.npz', nH_col_data=masked_nH_col.data, nH_col_mask=masked_nH_col.mask, mean_gz_data=masked_mean_gz.data, mean_gz_mask=masked_mean_gz.mask)

            lo_z_masked_nH_col=np.ma.masked_where(~lo_z_mask,masked_nH_col)
            lo_z_masked_mean_gz=np.ma.masked_where(~lo_z_mask,masked_mean_gz)

            np.savez(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked/lo_z_{mask}.npz', nH_col_data=lo_z_masked_nH_col.data, nH_col_mask=lo_z_masked_nH_col.mask, mean_gz_data=lo_z_masked_mean_gz.data, mean_gz_mask=lo_z_masked_mean_gz.mask)
    else:
        print('Masked Files Located')
    if 'plot_scatter' in kwargs: #Detects if plots are enabled
        if kwargs['plot_scatter']==True:
            if os.path.isdir(f'figures/{halo}/{bin_num}/nH_col_gz_scatter/{snap_num}')!=True:
                os.makedirs(f'figures/{halo}/{bin_num}/nH_col_gz_scatter/{snap_num}',exist_ok=True)
            plot.nH_col_gz_scatter(halo,bin_num,'xy',snap_num=snap_num) #Generates gas density plot
    elif 'all_plots' in kwargs: #Checks if override for all plots is enabled
        if kwargs['all_plots']==True:
            if os.path.isdir(f'figures/{halo}/{bin_num}/nH_col_gz_scatter/{snap_num}')!=True:
                os.makedirs(f'figures/{halo}/{bin_num}/nH_col_gz_scatter/{snap_num}',exist_ok=True)
            plot.nH_col_gz_scatter(halo,bin_num,'xy',snap_num=snap_num)


def get_threshold_behaviour_data(halo, **kwargs):
    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}') #Create halo directory if required
        print(f'Created directory for {halo}')
    else:
        print(f'\nDirectory for {halo} located')  

    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/threshold_behaviour')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/threshold_behaviour') #Create halo directory if required
        print(f'Created threshold vs time directory for {halo}')
    else:
        print(f'\nDirectory for {halo} located')  

    if 'snap_nums' not in kwargs:
        if 'redshifts' in kwargs:
            redshift_range=kwargs['redshifts']
            start_snap_num,start_snap_redshift=get_snap_num(halo,redshift_range[0])
            end_snap_num,end_snap_redshift=get_snap_num(halo,redshift_range[1])
            snap_nums=[start_snap_num,end_snap_num]
            snap_redshifts=[start_snap_redshift,end_snap_redshift]
        else:
            sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    else:
        snap_nums=kwargs['snap_nums']
        start_snap_redshift=get_redshift(halo,snap_nums[0])
        end_snap_redshift=get_redshift(halo,snap_nums[1])
        snap_redshifts=[start_snap_redshift,end_snap_redshift]

    for snap_num in snap_nums:
        snap_num=str(snap_num)
        while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
            snap_index=snap_nums.index(snap_num)
            snap_num='0'+snap_num
            snap_nums[snap_index]=snap_num

    if 'bins' in kwargs:
        bin_num=kwargs['bins']
    else:
        bin_num=512 

    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/threshold_behaviour/{bin_num}px')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/threshold_behaviour/{bin_num}px') #Create halo directory if required
        print(f'Created directory for {halo}')
    else:
        print(f'\nDirectory for {halo} located')  

    if os.path.isdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/threshold_behaviour/{bin_num}px/mass_fracs')!=True:
        os.mkdir(f'/cosma/apps/durham/dc-coll7/halos/{halo}/threshold_behaviour/{bin_num}px/mass_fracs') #Create halo directory if required


    snapshots=np.arange(int(snap_nums[0]),int(snap_nums[1])+1)

    all_data={'redshifts':[],'DLA':{'xy':[],'xz':[],'yz':[]},'lo_z_DLA':{'xy':[],'xz':[],'yz':[]},'subDLA':{'xy':[],'xz':[],'yz':[]},'lo_z_subDLA':{'xy':[],'xz':[],'yz':[]},'LymanLimit':{'xy':[],'xz':[],'yz':[]},'lo_z_LymanLimit':{'xy':[],'xz':[],'yz':[]}}

    for snap_num in snapshots:

        snap_num=str(snap_num)
        while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
            snap_num='0'+snap_num
        redshift=read_subfind_params(halo,snap_num=snap_num)['redshift']
        all_data['redshifts'].append(redshift)
        nH_col=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/nH_col.npy')
        mean_gz=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/mean_gz.npy')
        bin_masses=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/bin_masses.npy')
        
        planes=['xy','xz','yz']

        for plane in planes:
            plane_index=planes.index(plane)
            
            nH_col_proj=nH_col[plane_index]
            mean_gz_proj=mean_gz[plane_index]
            bin_masses_proj=bin_masses[plane_index]

            masks={'DLA':np.log10(nH_col_proj) >20.3,'subDLA': (20.3>np.log10(nH_col_proj)) & (np.log10(nH_col_proj) >19),'LymanLimit': (19>np.log10(nH_col_proj)) & (np.log10(nH_col_proj) >17.2)}

            total_mass=np.sum(bin_masses_proj.flatten())

            lo_z_mask=np.log10(mean_gz_proj)<-3

            for mask in masks:
                masked_bin_masses=np.ma.masked_where(~masks[mask],bin_masses_proj)
                lo_z_masked_bin_masses=np.ma.masked_where(~lo_z_mask,masked_bin_masses)

                masked_mass=np.sum(masked_bin_masses.compressed())
                lo_z_masked_mass=np.sum(lo_z_masked_bin_masses.compressed())

                masked_mass_frac=masked_mass/total_mass
                lo_z_masked_mass_frac=lo_z_masked_mass/total_mass

                all_data[mask][plane].append(masked_mass_frac)
                all_data[f'lo_z_{mask}'][plane].append(lo_z_masked_mass_frac)
    
    all_masks=['DLA','subDLA','LymanLimit','lo_z_DLA','lo_z_subDLA','lo_z_LymanLimit']

    for mask in all_masks:
        mask_data=np.array([all_data[mask]['xy'],all_data[mask]['xz'],all_data[mask]['yz']])
        np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/threshold_behaviour/{bin_num}px/mass_fracs/{mask}.npy',mask_data)
        
    np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/threshold_behaviour/redshifts.npy',all_data['redshifts'])
    

def get_stellar_masses(halo, **kwargs):
    if os.path.isdir(f'halos/{halo}')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'halos/{halo}') #Create halo directory if required
        print(f'Created directory for {halo}')
    else:
        print(f'\nDirectory for {halo} located')  

    if os.path.isdir(f'halos/{halo}/stellar_masses')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'halos/{halo}/stellar_masses') #Create halo directory if required
        print(f'Created stellar mass directory for {halo}')
    else:
        print(f'\nDirectory for {halo} located')

    if 'snap_nums' not in kwargs:
        if 'redshifts' in kwargs:
            redshift_range=kwargs['redshifts']
            start_snap_num,start_snap_redshift=get_snap_num(halo,redshift_range[0])
            end_snap_num,end_snap_redshift=get_snap_num(halo,redshift_range[1])
            snap_nums=[start_snap_num,end_snap_num]
            snap_redshifts=[start_snap_redshift,end_snap_redshift]
        else:
            sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    else:
        snap_nums=kwargs['snap_nums']
        start_snap_redshift=get_redshift(halo,snap_nums[0])
        end_snap_redshift=get_redshift(halo,snap_nums[1])
        snap_redshifts=[start_snap_redshift,end_snap_redshift]

    for snap_num in snap_nums:
        snap_num=str(snap_num)
        while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
            snap_index=snap_nums.index(snap_num)
            snap_num='0'+snap_num
            snap_nums[snap_index]=snap_num

    snapshots=np.arange(int(snap_nums[0]),int(snap_nums[1])+1)

    total_stellar_masses=[]
    redshifts=[]

    for snap_num in snapshots:
        snap_num=str(snap_num)
        while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
            snap_num='0'+snap_num

        stellar_masses=read_raw_file(halo,'stars','mass',snap_num=snap_num)
        total_stellar_mass=np.sum(stellar_masses)
        total_stellar_masses.append(total_stellar_mass.value)

        redshift=get_redshift(halo,snap_num)
        redshifts.append(redshift)
    
    total_stellar_masses=np.array(total_stellar_masses)
    redshifts=np.array(redshifts)

    np.save(f'halos/{halo}/stellar_masses/redshifts.npy',redshifts)
    np.save(f'halos/{halo}/stellar_masses/stellar_masses.npy',total_stellar_masses)

def get_properties_table(halo):
    loc='/cosma8/data/dp004/lyra/original_sample/' #Uses default location if none provided
    suffix='/output/' #Suffix to ensure data from each snapshot is read
    
    final_snapshot=get_snap_nums(halo)[-1]

    subfind_data = ar.gadget_subfind.load_subfind(final_snapshot, dir=loc + halo + suffix)
       

    
  
                

                

              
    
        



    
        










    

        
    
    
            
                
    
    


            



    

    
    

    
            
            
    
        


            


           




                

    