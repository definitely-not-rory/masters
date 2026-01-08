#---Imports and Admin---
from imports import *
from sim_folder_readers import *

def create_redshift_file(halo): #Function to create snapshot-redshift files for given simulations
    if os.path.exists('sim_snapshot_redshifts/'+halo+'.npy'): #Check if redshift file has already been created
            print(halo+' Redshift File Already Exists')
    else: #For simulations without created redshift files
        loc = '/cosma8/data/dp004/lyra/original_sample/' #Assign target data storage directory
        suffix = '/output/'
        
        snap_nums = get_sn_nums(halo) #Retrieve available snapshots
        
        def get_redshift(snap_num):
            subfind = ar.gadget_subfind.load_subfind(int(snap_num), dir=loc + halo + suffix, onlyHeader=True) #Import subfind header only for given snapshot number
            redshift=subfind.redshift #Retrieve snapshot redshift value from header
            return redshift
        
        snap_redshifts=[get_redshift(snap_num) for snap_num in snap_nums] #Pull all redshifts for available snapshots in chosen simulation directory
        
        storage_array=np.array([snap_nums,snap_redshifts]) #Create 2D storage array for saving to external file
        np.save('sim_snapshot_redshifts/'+halo+'.npy',storage_array) #Save to external file

def get_redshift_snap_num(halo,redshift):#Function to return the closest available snapshot number for a given halo and redshift
    data=np.load('sim_snapshot_redshifts/'+halo+'.npy') #Import given halos data
    index=np.argmin(np.abs(np.float64(data[1])-redshift)) #Locate index of redshift array entry with smallest absolute difference to target redshift
    target_snap_num=data[0,index] #Retrieve correct snapshot number from array
    return target_snap_num