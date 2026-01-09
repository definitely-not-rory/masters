from imports import *
from halo_readers import get_halos, get_redshifts, get_snap_num, get_redshift
from npy_data_readers import read_raw_file, read_subfind_params
import plot_generation as plot 
import processing as calc  

#test

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
    
    if os.path.isdir(f'halos/{halo}')!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(f'halos/{halo}') #Create halo directory if required
        print(f'Created directory for {halo}')
    else:
        print(f'\nDirectory for {halo} located')

    #--- Snapshot Number Selection ---

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
            print(f'\nTarget Redshift: 0, Using Final Snapshot ({snap_num})')
        else:
            snap_num,snap_redshift=get_snap_num(halo,target_redshift)
    else:
        sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")') #Error message for incorrect snapshot number/redshift

    if 'plot_dz_snap' in kwargs: #Detects if plots are enabled
        if kwargs['plot_dz_snap']==True:
            if 'redshift' in kwargs:
                plot.dz_snapshot(halo,snap_num,target_redshift=target_redshift,snap_redshift=snap_redshift) #Generates dz vs snapshot plot with target redshift elements if redshift was used instead of snapshot number
            else:
                plot.dz_snapshot(halo,snap_num) #Generates dz vs snapshot plot without target redshift elements if snapshot number was used
    elif 'all_plots' in kwargs:
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

    if os.path.isdir(f'halos/{halo}/{snap_num}')!=True: #Checks for dedicated directory for chosen snapshot
        print(f'\nNo Directory Detected for {halo} Snapshot {snap_num}') 
        os.mkdir(f'halos/{halo}/{snap_num}') #Creates snapshot directory if necessary
        print(f'{halo} Snapshot {snap_num} Directory Created')
    else:
        print(f'\n{halo} Snapshot {snap_num} Directory Located')    

    #--- Subfind Data Handling ---

    if os.path.isdir(f'halos/{halo}/{snap_num}/subfind')!=True: #Checks for subfind data directory within snapshot directory
        print(f'\nNo {halo} Snapshot {snap_num} subfind directory detected')
        os.mkdir(f'halos/{halo}/{snap_num}/subfind') #Creates subfind directory if necessary
        print(f'Subfind directory for {halo} Snapshot {snap_num} created')
    else:
        print(f'\n{halo} Snapshot {snap_num} Subfind Directory Located')

    print(f'\nImporting {halo} Snapshot {snap_num} Subfind Data')
    subfind_data = ar.gadget_subfind.load_subfind(int(snap_num), dir=loc + halo + suffix) #Import subfind dataset
    print('Subfind Data imported')

    if os.path.exists(f'halos/{halo}/{snap_num}/subfind/fof_positions.npy')!=True or os.path.exists(f'halos/{halo}/{snap_num}/subfind/fof_masses.npy')!=True or os.path.exists(f'halos/{halo}/{snap_num}/subfind/halo_params.npy')!=True: #Detects if subfind data arrays have been already imported and stored external
        print(f'\nNo Subfind FoF Data Located\Loading subfind data for {halo} from Snapshot {snap_num}') #Creates relevant external data arrays if required

        #--- Subfind-Wide Data Processing ---
        sf_positions=subfind_data.data['fpos']  #FoF Positions (in Mpc)
        all_sf_masses=subfind_data.data['fmty'] #FoF Masses Across All Types (in 10^10 M_sol)

        high_res_mask=(all_sf_masses[:,2]==0) & (all_sf_masses[:,3]==0) #Masks mass data to exclude low-res structures
        high_res_sf_masses=all_sf_masses[high_res_mask]
        sf_masses=np.sum(high_res_sf_masses,axis=1) #Calculates total high-res mass for each FoF group
        
        np.save(f'halos/{halo}/{snap_num}/subfind/fof_positions.npy',sf_positions) #Saves subfind-wide arrays to external .npy files
        print('Subfind FoF position data saved')

        np.save(f'halos/{halo}/{snap_num}/subfind/fof_masses.npy',sf_masses)
        print('Subfind FoF mass data saved')

        #--- FoF Group 1/Main Halo Data Processing ---
        halo_pos=sf_positions[0] #Main Halo Position (in Mpc)
        halo_mass=sf_masses[0] #Main Halo Mass (in 10^10 M_sol)
        halo_r200=subfind_data.data['frc2'][0] #Main Halo R_200_crit (in Mpc)

        redshift=subfind_data.redshift #Retrieve snapshot redshift value from header
        
        saved_params=np.array([redshift,halo_pos,halo_mass,halo_r200],dtype='object') #Create storage array for FoF Group 1 parameters
        
        np.save(f'halos/{halo}/{snap_num}/subfind/halo_params.npy',saved_params) #Save FoF Group 1 parameters to external .npy file
        print('Subfind Halo Parameters Saved')
    else:
        print('\nFoF Subfind Data Located')
    
    if 'plot_fof_scatter' in kwargs: #Detects if plots are enabled
        if kwargs['plot_fof_scatter']==True:
            plot.fof_scatter(halo,snap_num=snap_num) #Generates FoF scatter plot
    elif 'all_plots' in kwargs: #Checks if override for all plots is enabled
        if kwargs['all_plots']==True:
            plot.fof_scatter(halo,snap_num=snap_num)

    #--- Snapshot Data Handling ---

    #--- Directory Management ---
    if os.path.isdir(f'halos/{halo}/{snap_num}/raw')!=True: #Checks for raw snapshot data directory within snapshot directory
        print(f'\nNo {halo} raw Snapshot {snap_num} data directory detected')
        os.mkdir(f'halos/{halo}/{snap_num}/raw') #Creates raw data directory if necessary

        os.mkdir(f'halos/{halo}/{snap_num}/raw/gas') #Creates raw matter type directories
        os.mkdir(f'halos/{halo}/{snap_num}/raw/dm')
        os.mkdir(f'halos/{halo}/{snap_num}/raw/stars')

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

    current_gas_data=[file[:-4] for file in os.listdir(f'halos/{halo}/{snap_num}/raw/gas')] #Import list of current files saved for each file
    current_dm_data=[file[:-4] for file in os.listdir(f'halos/{halo}/{snap_num}/raw/dm')]
    current_stars_data=[file[:-4] for file in os.listdir(f'halos/{halo}/{snap_num}/raw/stars')]

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

                np.save(f'halos/{halo}/{snap_num}/raw/dm/pos.npy',dm_pos) #Saves DM position data to external .npy file
                print(f'DM pos data saved externally')
                
                #Loads necessary parameters to determine constant DM particle mass
                h=snapshot_data.header['HubbleParam'] #Import h and header DM particle mass
                dm_particle_mass=snapshot_data.header['MassTable'][1]/h #Apply relevant conversion from header DM particle mass to 'real' mass
                print(f'\nDM Mass Parameters data loaded')

                dm_params=np.array([h,dm_particle_mass]) #Stores all DM mass parameters in array for external file
                np.save(f'halos/{halo}/{snap_num}/raw/dm/dm_params.npy',dm_params) #Saves DM parameters to external .npy file
                print(f'DM Mass Parameters data saved externally')
                
                dm_mass=np.full_like(dm_pos[:,0],dm_particle_mass) #Creates 'normal' DM mass array using params and loaded snapshot data
                print(f'\nDM mass data calculated')

                np.save(f'halos/{halo}/{snap_num}/raw/dm/mass.npy',dm_mass) #Saves DM mass array
                print(f'DM mass data saved externally')     
            else: #Completes all non-exception exports
                to_save=req_snapshot_data[type_index] #Loads all requires snapshot data labels
                for param in to_save: #Exports each required parameter for target type
                    data=snapshot_data.data[param] #Loads target type's target parameter from snapshot
                    print(f'\n{type_name} {param} data loaded')

                    directory=type_directories[type_index] #Finds correct type directory for data storage
                    np.save(f'halos/{halo}/{snap_num}/raw/{directory}/{param}.npy',data) #Saves parameter data to external .npy file
                    print(f'{type_name} {param} data saved externally')
        else:
            print(f'\n{halo} Snapshot {snap_num} {type_name} raw data complete')
    print(f'\n{halo} Snapshot {snap_num} Raw Data Import Complete')

def get_mass_density_data(halo,**kwargs):
    if 'snap_num' not in kwargs:
        if 'redshift' in kwargs:
            target_redshift=kwargs['redshift']
            snap_num,snap_redshift=get_snap_num(halo,target_redshift)
        else:
            sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    else:
        snap_num=kwargs['snap_num']
        snap_redshift=get_redshift(halo,snap_num)

    halo_pos=read_subfind_params(halo,snap_num=snap_num)['halo_pos']
    redshift=read_subfind_params(halo,snap_num=snap_num)['redshift']

    type_directories=['gas','dm','stars']
    params=[['pos','mass'],['pos','mass'],['pos','mass']]
    type_names=['Gas','DM','Stars']
    
    loaded_data = {}

    for matter_type in type_directories:
        type_index=type_directories.index(matter_type)
        type_name=type_names[type_index]
        to_load=params[type_index]

        if os.path.exists(f'halos/{halo}/{snap_num}/raw/{matter_type}/pos.npy')!=True or os.path.exists(f'halos/{halo}/{snap_num}/raw/{matter_type}/mass.npy')!=True:
            print(f'\nSnapshot {snap_num} {type_name} Raw Data Incomplete')
            get_raw_data(halo,snap_num=snap_num)
        else:
            print(f'\nSnapshot {snap_num} {type_name} raw data located')
        
        type_data={param:read_raw_file(halo,matter_type,param,snap_num=snap_num) for param in to_load}
        
        loaded_data[matter_type]=type_data

    for matter_type, data in loaded_data.items():
        type_name=type_names[type_directories.index(matter_type)]

        if os.path.exists(f'halos/{halo}/{snap_num}/raw/{matter_type}/rel_pos.npy')!=True:
            print(f'{type_name} Relative Position Data Not Found, Generating')

            rel_pos=calc.to_rel(data['pos'],halo_pos)
            loaded_data[matter_type]['rel_pos']=rel_pos

            np.save(f'halos/{halo}/{snap_num}/raw/{matter_type}/rel_pos.npy',rel_pos.value)
            print(f'{type_name} Relative Position File Generated')
        else:
            rel_pos=read_raw_file(halo,matter_type,'rel_pos',snap_num=snap_num)
            print(f'{type_name} Relative Positions Loaded')

            loaded_data[matter_type]['rel_pos']=rel_pos

    for matter_type, data in loaded_data.items():
        type_name=type_names[type_directories.index(matter_type)]
        if os.path.exists(f'halos/{halo}/{snap_num}/raw/{matter_type}/radii.npy')!=True:
            print(f'{type_name} Radial Position Data Not Found, Generating')

            radii=calc.to_rad(data['pos'],halo_pos)
            loaded_data[matter_type]['radii']=radii

            np.save(f'halos/{halo}/{snap_num}/raw/{matter_type}/radii.npy',radii.value)
            print(f'{type_name} Radial Position File Generated')
        else:
            radii=read_raw_file(halo,matter_type,'radii',snap_num=snap_num)
            print(f'{type_name} Radial Positions Loaded')

            loaded_data[matter_type]['radii']=radii
    
    if os.path.isdir(f'halos/{halo}/{snap_num}/binned')!=True:
        os.mkdir(f'halos/{halo}/{snap_num}/binned')
                      
    if 'bins' in kwargs:
        bin_num=kwargs['bins']
    else:
        bin_num=512

    if os.path.isdir(f'halos/{halo}/{snap_num}/binned/{bin_num}')!=True:
        os.mkdir(f'halos/{halo}/{snap_num}/binned/{bin_num}')
    
    if os.path.isdir(f'halos/{halo}/{snap_num}/binned/{bin_num}/total_mass')!=True:
        os.mkdir(f'halos/{halo}/{snap_num}/binned/{bin_num}/total_mass')
    
    for matter_type, data in loaded_data.items():
        type_name=type_names[type_directories.index(matter_type)]
        if os.path.exists(f'halos/{halo}/{snap_num}/binned/{bin_num}/total_mass/{matter_type}.npy')!=True:
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
            np.save(f'halos/{halo}/{snap_num}/binned/{bin_num}/total_mass/{matter_type}.npy',save_data)
            print(f'\n{matter_type} projected density file created')
        else:
            loaded_data[matter_type]['proj_dens']=np.load(f'halos/{halo}/{snap_num}/binned/{bin_num}/total_mass/{matter_type}.npy')
            print(f'\n{matter_type} projected density file loaded')
            








'''xy_binned_mass=stats.binned_statistic_2d(pos[:,0].to_value(units.cm),pos[:,1].to_value(units.cm),mass.to_value(units.g)/factor,bins=[bin_num,bin_num],statistic='sum')
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

    return xy_binned_mass.statistic/binarea_xy,xz_binned_mass.statistic/binarea_xz,yz_binned_mass.statistic/binarea_yz'''

    
    


            



    

    
    

    
            
            
    
        


            


           




                

    