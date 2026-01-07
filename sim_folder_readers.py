from imports import *

def get_simulation_list(): #Function to return all simulation directories
    folders = [folder for folder in os.listdir('/cosma8/data/dp004/lyra/original_sample/') if os.path.isdir(os.path.join('/cosma8/data/dp004/lyra/original_sample/', folder))] #Retrieve list of directories
    sim_folders = [folder for folder in folders if 'DMO' not in folder] #Select only simulation directories
    output_str = '\n'.join(run for run in sim_folders) #Append all simulation directories to a string
    return output_str

def find_last_snapshot(name): #Function to return last snapshot in a given simulation directory
    folders = [folder for folder in os.listdir('/cosma8/data/dp004/lyra/original_sample/'+name+'/output') if os.path.isdir(os.path.join('/cosma8/data/dp004/lyra/original_sample/'+name+'/output', folder))] #Return all directories in given simulation directory
    group_folders = [folder for folder in folders if folder.startswith('groups_')] #Select only directories with snapshot data
    group_numbers=[int(name[7:]) for name in group_folders] #Return numbers from ends of snapshot data directories
    last_group=(max(group_numbers)) #Select highest directory number
    return last_group

def get_sn_nums(name): #Function to return all available snapshot numbers for a given simulation directory
    folders = [folder for folder in os.listdir('/cosma8/data/dp004/lyra/original_sample/'+name+'/output') if os.path.isdir(os.path.join('/cosma8/data/dp004/lyra/original_sample/'+name+'/output', folder))] #Return all directories in given simulation directory
    group_folders = [folder for folder in folders if folder.startswith('groups_')] #Select only directories with snapshot data
    group_numbers=[name[7:] for name in group_folders] #Return numbers from ends of snapshot data directories
    return sorted(group_numbers)[5:]