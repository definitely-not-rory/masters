from imports import *

def get_halos(**alt_path): #Function to return all halo directory names
    if len(alt_path.items())!=0: #Detects if alternate directory for halo data is provided
        path=alt_path['path']
        print(f'Alternative Directory Selected: {path}')
    else:
        path='/cosma8/data/dp004/lyra/original_sample/' #Selects default directory if not other path is given
        print(f'Using default directory: {path}')
    folders = [folder for folder in os.listdir(path) if os.path.isdir(os.path.join(path, folder))] #Retrieve list of directories
    sim_folders = [folder for folder in folders if 'DMO' not in folder] #Selects only non-Dark-Matter-Only directories
    return sim_folders

def get_snap_nums(halo,**alt_path): #Function to return all available snapshot numbers for a given halo
    if len(alt_path.items())!=0: #Check for alternate directory
        path=alt_path['path']
    else:
        path='/cosma8/data/dp004/lyra/original_sample/' #Assume default data location/directory
        
    folders = [folder for folder in os.listdir(f'{path}{halo}/output') if os.path.isdir(os.path.join(f'{path}{halo}/output', folder))] #Return all directories in given halo's directory
    group_folders = [folder for folder in folders if folder.startswith('groups_')] #Select only directories with snapshot data
    group_numbers=sorted([name[7:] for name in group_folders]) #Return ordered numbers from ends of snapshot data directories
    return group_numbers

def get_redshift(halo,snap_num,**kwargs): #Auxilliary function to load one individual redshift for a given snapshot
    snap_num=str(snap_num) #Loads snapshot number integer input
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
        snap_num='0'+snap_num
    if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/redshifts.npy')==True: #Detects if redshift data file already exists for provided halo
        print(f'\n Redshift file for {halo} located')
        snap_nums,snap_redshifts=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/redshifts.npy') #Import existing redshift-snapshot data
        snap_index=snap_nums.tolist().index(snap_num) #Find index of requested snapshot
        redshift=np.float64(snap_redshifts[snap_index]) #Select redshift of desired snapshot
        print(f'\nSnapshot Number: {snap_num}\nSnapshot Redshift: {redshift}')
        return redshift
    else:
        if 'loc' not in kwargs:
            loc='/cosma8/data/dp004/lyra/original_sample/'
            print(f'Using default directory: {loc}')
        else:
            loc = kwargs['loc'] #Selects requested halo directory
        suffix='/output/' #Suffix to ensure data from each snapshot is read
        subfind = ar.gadget_subfind.load_subfind(int(snap_num), dir=loc + halo + suffix, onlyHeader=True) #Import subfind header only for given snapshot number
        redshift=subfind.redshift #Retrieve snapshot redshift value from header
        if 'print_data' in kwargs: #Detects if a printed output is requested
            print_text=kwargs['print_data']
            if print_text==True:
                print(f'\nSnapshot Number: {snap_num}\nSnapshot Redshift: {redshift}') #Prints output if necessary
        return redshift

def get_redshifts(halo,**kwargs): #Function to import and store the redshift data of all snapshots for a given halo
    if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/redshifts.npy')!=True: #Detects if redshift data file already exists for provided halo
        print(f'No redshift file located for {halo}, generating')
        if 'path' in kwargs: #Detects if alternate halo data directory has been provided
            loc=kwargs['path']
        else:
            loc='/cosma8/data/dp004/lyra/original_sample/' #Uses default location if none provided
        
        snap_nums=get_snap_nums(halo) #Get all available snapshot numbers from halo's directory
        
        snap_redshifts=[get_redshift(halo,snap_num,loc=loc) for snap_num in snap_nums] #Pull all redshifts for available snapshots in chosen simulation directory
        save_data=np.array([snap_nums,snap_redshifts]) #Create 2D storage array for saving to external file
        np.save(f'/cosma/apps/durham/dc-coll7/halos/{halo}/redshifts.npy',save_data) #Saves data to .npy file in halo's .npy storage directory
    else: #If redshift file already exists
        print(f'{halo} Redshift File Located') 
        snap_nums,snap_redshifts=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/redshifts.npy') #Load .npy stored data
    return snap_nums, snap_redshifts

def get_snap_num(halo,target_redshift):
    if os.path.exists(f'/cosma/apps/durham/dc-coll7/halos/{halo}/redshifts.npy')!=True: #Detects if redshift data file already exists for provided halo
        print(f'\nNo redshift file for {halo}, generating...')
        snap_nums,snap_redshifts=get_redshifts(halo)
    else:
        snap_nums,snap_redshifts=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/redshifts.npy')
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

