from imports import *
from halo_readers import get_halos, get_snap_nums, get_redshifts   
import plot_generation as plot   

def get_raw_data(halo,**kwargs): #Function to import all raw subfind and snapshot data into .npy files for a given halo and snapshot number/redshift
    if 'path' in kwargs: #Detects for provision of alternate data directory
        alt_path=kwargs['path']
        available_halos=get_halos(path=alt_path) #Loads halos in alternate directory
    else:
        available_halos=get_halos() #Loads halos in default directory if no alternate is provided
    
    if halo not in available_halos: #Detects if requested halo exists
        print('Halo does not exist, please select from:') #Returns lis of available halo in provided/default directory
        print(*available_halos,sep='\n')
        sys.exit('Halo Name Not Found in Directory') #Error message for incorrect halo name
    
    if os.path.isdir(halo)!=True: #Detect for dedicated .npy storage directory for halo raw data
        print('\nHalo Directory Not Present')
        os.mkdir(halo) #Create halo directory if required
        print(f'Created directory for {halo}')
    else:
        print(f'\nDirectory for {halo} located')

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
            index=np.argmin(np.abs(np.float64(snap_redshifts)-target_redshift)) #Locate index of redshift array entry with smallest absolute difference to target redshift
            snap_num=snap_nums[index] #Retrieve correct snapshot number from array
            print(f'\n--- Target Redshift: {target_redshift} ---\nSnapshot Number Selected: {snap_num}\nSnapshot True Redshift: {snap_redshifts[index]}\nRedshift Difference: {np.abs(snap_redshifts[index]-target_redshift)/target_redshift*100}%') #Print % error on target vs snapshot redshift       
    else:
        sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")') #Error message for incorrect snapshot number/redshift
  
    if 'plot_dz_snap' in kwargs: #Detects if plots are enabled
        if kwargs['plot_dz_snap']==True:
            if 'redshift' in kwargs:
                plot.dz_snapshot(halo,snap_num,target_redshift=target_redshift,snap_redshift=snap_redshifts[index]) #Generates dz vs snapshot plot with target redshift elements if redshift was used instead of snapshot number
            else:
                plot.dz_snapshot(halo,snap_num) #Generates dz vs snapshot plot without target redshift elements if snapshot number was used

    if 'path' in kwargs: #Detects for provision of alternate data directory
        loc=kwargs['path']
    else:
        loc='/cosma8/data/dp004/lyra/original_sample/' #Uses default location if none provided
    suffix='/output/' #Suffix to ensure data from each snapshot is read

    if os.path.isdir(f'{halo}/{snap_num}')!=True:
        print(f'\nNo Directory Detected for {halo} Snapshot {snap_num}')
        os.mkdir(f'{halo}/{snap_num}')
        print(f'{halo} Snapshot {snap_num} Directory Created')
    else:
        print(f'\n{halo} Snapshot {snap_num} Directory Located')    


    if os.path.isdir(f'{halo}/{snap_num}/subfind')!=True:
        print(f'\nNo {halo} Snapshot {snap_num} subfind directory detected')
        os.mkdir(f'{halo}/{snap_num}/subfind')
        print(f'Subfind directory for {halo} Snapshot {snap_num} created')
    else:
        print(f'\n{halo} Snapshot {snap_num} Subfind Directory Located')

    if os.path.exists(f'{halo}/{snap_num}/subfind/fof_positions.npy')!=True or os.path.exists(f'{halo}/{snap_num}/subfind/fof_masses.npy')!=True or os.path.exists(f'{halo}/{snap_num}/subfind/halo_params.npy')!=True:
        print(f'\nNo Subfind FoF Positions Data Located\nImporting subfind data for {halo} from Snapshot {snap_num}')
        subfind_data = ar.gadget_subfind.load_subfind(int(snap_num), dir=loc + halo + suffix) #Import subfind dataset
        print('Subfind data imported')

        sf_positions=subfind_data.data['fpos'] #Imports necessary subfind data
        all_sf_masses=subfind_data.data['fmty']

        high_res_mask=(all_sf_masses[:,2]==0) & (all_sf_masses[:,3]==0)
        high_res_sf_masses=all_sf_masses[high_res_mask]
        sf_masses=np.sum(high_res_sf_masses,axis=1)
        
        np.save(f'{halo}/{snap_num}/subfind/fof_positions.npy',sf_positions)
        print('Subfind FoF position data saved')

        np.save(f'{halo}/{snap_num}/subfind/fof_masses.npy',sf_masses)
        print('Subfind FoF mass data saved')

        halo_pos=sf_positions[0]
        halo_mass=sf_masses[0]
        halo_r200=subfind_data.data['frc2'][0]
        redshift=subfind_data.redshift #Retrieve snapshot redshift value from header
        
        saved_params=np.array([redshift,halo_pos,halo_mass,halo_r200],dtype='object')
        
        np.save(f'{halo}/{snap_num}/subfind/halo_params.npy',saved_params)
        print('Subfind Halo Parameters Saved')
    else:
        print('\nSubfind Data Located')
    
    if 'plot_fof_scatter' in kwargs: #Detects if plots are enabled
        if kwargs['plot_fof_scatter']==True:
            plot.fof_scatter(halo,snap_num) #Generates FoF scatter plot


           





                

    