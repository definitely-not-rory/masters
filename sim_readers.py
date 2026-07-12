from imports import *

def get_all_snapshots(halo,path='/cosma8/data/dp004/lyra/original_sample/'): #Function to return all available snapshot numbers for a given halo
    if os.path.exists(os.path.join(f'{path}{halo}/output')):   #Ensures selected directory contains simulation snapshot data
        folders = [folder for folder in os.listdir(f'{path}{halo}/output') if os.path.isdir(os.path.join(f'{path}{halo}/output', folder))] #Return all directories in given halo's directory
        
        group_folders = [folder for folder in folders if folder.startswith('groups_')] #Select only directories with snapshot data
        group_numbers=sorted([name[7:] for name in group_folders]) #Return ordered numbers from ends of snapshot data directories
        
        if len(group_numbers)>1: #Ensures there is sufficient snapshots in selected halo dataset
            print(f'Located {len(group_numbers)} snapshots for {halo} in {path}.\n')
            return group_numbers
        else:
            print(f'Insufficient ({len(group_numbers)}) snapshots for {halo}.\n')
            return False
    else:
        print(f'{path}{halo} is not a LYRA simulation directory.\n')
        return False

def get_sims(path='/cosma8/data/dp004/lyra/original_sample/'): #Function to return all halo directory names
    
    if path!='/cosma8/data/dp004/lyra/original_sample/': #Detects if alternate directory for halo data is provided
        print(f'Alternative Directory Selected: {path}\n')
    else:
        print(f'Using default directory: {path}\n')
    
    folders = [folder for folder in os.listdir(path) if os.path.isdir(os.path.join(path, folder))] #Retrieve list of directories
    non_DMO_folders = [folder for folder in folders if 'DMO' not in folder] #Selects only non-Dark-Matter-Only directories
    
    sim_folders=[folder for folder in non_DMO_folders if get_all_snapshots(folder)!=False]

    print('Simulations With Sufficient Snapshots:')
    for sim in sim_folders:
        print(f'- {sim}')

    return sim_folders

def get_redshift(halo,snap_num,data_folder='/cosma8/data/dp004/lyra/original_sample/',results_folder='/cosma/apps/durham/dc-coll7/halos/'): #Auxilliary function to load one individual redshift for a given snapshot
    
    snap_num=str(snap_num) #Loads snapshot number integer input
    
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
        snap_num='0'+snap_num
    
    if os.path.exists(f'{results_folder}{halo}/redshifts.npy')==True: #Detects if redshift data file already exists for provided halo
        
        print(f'\n Existing redshifts.npy file for {halo} located in {results_folder}.')
        
        snap_nums,snap_redshifts=np.load(f'{results_folder}{halo}/redshifts.npy') #Import existing redshift-snapshot data
        snap_index=snap_nums.tolist().index(snap_num) #Find index of requested snapshot
        
        redshift=np.float64(snap_redshifts[snap_index]) #Select redshift of desired snapshot
        
        print(f'\nSnapshot Number: {snap_num}\nSnapshot Redshift: {redshift}')
        return redshift
    else:
        print(f'redshifts.npy not found in {results_folder}{halo}, loading redshift from Subfind.\n')

        if data_folder !='/cosma8/data/dp004/lyra/original_sample/':
            if os.path.exists(f'{data_folder}{halo}')==True: 
                print(f'Alternative Data Directory Selected: {data_folder}')
            else:
                print(f'No data for {halo} found at {data_folder}')
                return False
        else:
            print(f'Loading from default data directory: {data_folder}')
        
        suffix='/output/' #Suffix to ensure data from each snapshot is read
        
        subfind = ar.gadget_subfind.load_subfind(int(snap_num), dir=data_folder + halo + suffix, onlyHeader=True) #Import subfind header only for given snapshot number
        redshift=subfind.redshift #Retrieve snapshot redshift value from header
        
        print(f'\nSnapshot Number: {snap_num}\nSnapshot Redshift: {redshift}') #Prints output if necessary
        return redshift

def get_redshifts(halo,data_folder='/cosma8/data/dp004/lyra/original_sample/',results_folder='/cosma/apps/durham/dc-coll7/halos/'): #Function to import and store the redshift data of all snapshots for a given halo
    if os.path.exists(f'{results_folder}{halo}/redshifts.npy')!=True: #Detects if redshift data file already exists for provided halo
        print(f'No redshift file located for {halo} in {results_folder}, generating at {results_folder}{halo}/redshifts.npy.')
        
        snap_nums=get_all_snapshots(halo,path=data_folder) #Get all available snapshot numbers from halo's directory
        snap_redshifts=[get_redshift(halo,snap_num,data_folder=data_folder,results_folder=results_folder) for snap_num in snap_nums] #Pull all redshifts for available snapshots in chosen simulation directory
        save_data=np.array([snap_nums,snap_redshifts]) #Create 2D storage array for saving to external file

        if os.path.exists(f'{results_folder}{halo}')!=True:
            
            valid_input=False #Yes/No input auxilliary variables.
            valid_inputs=['y','Y','n','N']
        
            while valid_input==False: #Only progress function until Yes/No receives input of 'Y', 'y', 'N' or 'n'.
                overwrite=input(f'Halo directory for {halo} does not exist in {results_folder}, would you like to create one? (y/n):') 
                if overwrite in valid_inputs: #Break loop if a valid input is received.
                    valid_input=True
                    if overwrite in valid_inputs[2:]: #Check if overwrite permission is 'No', and halt function if it is.
                        return
                else:
                    print(f'\nCreating results directory for {halo} at {results_folder}')
        
        os.makedirs(f'{results_folder}{halo}',exist_ok=True)
        np.save(f'{results_folder}{halo}/redshifts.npy',save_data) #Saves data to .npy file in halo's .npy storage directory
        print(f'Snapshot redshift data for {halo} save at {results_folder}/redshifts.npy')
    else: #If redshift file already exists
        print(f'{halo} Redshift File Located in {results_folder}.') 
        snap_nums,snap_redshifts=np.load(f'{results_folder}{halo}/redshifts.npy') #Load .npy stored data
    
    return snap_nums, snap_redshifts

def get_snap_num(halo,target_redshift,data_folder='/cosma8/data/dp004/lyra/original_sample/',results_folder='/cosma/apps/durham/dc-coll7/halos/'):
    if os.path.exists(f'{results_folder}{halo}/redshifts.npy')!=True: #Detects if redshift data file already exists for provided halo
        print(f'\nNo redshift file for {halo} in {results_folder}, generating...')
        snap_nums,snap_redshifts=get_redshifts(halo,data_folder=data_folder,results_folder=results_folder)
    else:
        snap_nums,snap_redshifts=np.load(f'{results_folder}{halo}/redshifts.npy')
        print(f'\n {halo} redshift file imported')
    
    if target_redshift==0:
        snap_num=snap_nums[-1]
        snap_redshift=np.float64(snap_redshifts[-1])
    else:
        index=np.argmin(np.abs(np.float64(snap_redshifts)-np.float64(target_redshift))) #Locate index of redshift array entry with smallest absolute difference to target redshift
        snap_num=snap_nums[index] #Retrieve correct snapshot number from array
        snap_redshift=np.float64(snap_redshifts[index])
    
    print(f'\n--- Target Redshift: {target_redshift} ---\nSnapshot Number Selected: {snap_num}\nSnapshot True Redshift: {snap_redshift}\nRedshift Difference: {np.abs(snap_redshift-np.float64(target_redshift))/np.float64(target_redshift)*100}%') #Print % error on target vs snapshot redshift
    return snap_num, snap_redshift

