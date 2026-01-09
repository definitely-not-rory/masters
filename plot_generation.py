from imports import *
from halo_readers import get_redshift, get_snap_num
from npy_data_readers import read_raw_file, read_subfind_params
import processing as calc

def dz_snapshot(halo,snap_num,**redshifts):
    snap_nums,snap_redshifts=np.float64(np.load(f'halos/{halo}/redshifts.npy'))

    dz=np.abs(np.diff(snap_redshifts))

    fig_dz_snapshot,ax_dz_snapshot=plt.subplots()

    ax_dz_snapshot.plot(snap_nums[1:],dz,c='r')
    ax_dz_snapshot.axvline(int(snap_num),c='k',linestyle='dashed')

    ax_dz_snapshot.set_xlabel('Snapshot Number')
    ax_dz_snapshot.set_ylabel('Absolute Change In Redshift ($\Delta z$)')

    plt.text(snap_nums[1]+2, dz[0]-.2, '$\Delta z$', fontsize=14,c='r')

    z_snap_num=ax_dz_snapshot.twinx()
    z_snap_num.plot(snap_nums,snap_redshifts,c='b',alpha=.5)
    z_snap_num.set_ylabel('Snapshot Redshift ($z$)')

    plt.text(snap_nums[0]-5, snap_redshifts[0]+.3, '$z$', fontsize=14,c='b')
    plt.text(int(snap_num)+2, 40, f'Snapshot Number $= {snap_num}$', fontsize=12,c='k',rotation='vertical')

    display_redshift=np.round(get_redshift(halo,snap_num),3)
    display_halo=halo.replace('_',' ')

    plt.title(f'{display_halo}, $z = {display_redshift}$')

    if 'target_redshift' in redshifts:
        target_redshift=redshifts['target_redshift']
        snap_redshift=redshifts['snap_redshift']

        z_snap_num.axhline(target_redshift,c='b',linestyle='dashed')
        z_snap_num.axhline(snap_redshift,c='b')

        zoomed=ax_dz_snapshot.inset_axes([0.2,0.35,0.45,0.5])
        zoomed.plot(snap_nums[1:],dz,c='r')
        zoomed.axvline(int(snap_num),c='k',linestyle='dashed')

        zoomed.set_xlim(int(snap_num)-10,int(snap_num)+10)

        zoomed.set_xlabel('Snapshot Number',backgroundcolor='white')
        zoomed.set_ylabel('Absolute Change In Redshift ($\Delta z$)')

        zoomed_z_snap_num=zoomed.twinx()
        zoomed_z_snap_num.plot(snap_nums,snap_redshifts,c='b',alpha=.5)
        zoomed_z_snap_num.axhline(target_redshift,c='b',linestyle='dashed')
        zoomed_z_snap_num.axhline(snap_redshift,c='k')

        plt.text(int(snap_num)-3, target_redshift-.4, f'Snapshot Number $= {snap_num}$', fontsize=10,c='k',rotation='vertical',backgroundcolor='white')
        plt.text(int(snap_num)+5,target_redshift+.05,f'Target\n Redshift $= {target_redshift}$',c='b',ha='center')
        plt.text(int(snap_num)+5,snap_redshift-.2,f'Snapshot\n Redshift $= {np.round(snap_redshift,3)}$',ha='center')


        zoomed_z_snap_num.set_ylim(target_redshift-.5,target_redshift+.5)

        zoomed_z_snap_num.set_ylabel('Snapshot Redshift ($z$)',backgroundcolor='white')

        z_snap_num.indicate_inset_zoom(zoomed,edgecolor="black",alpha=1)

    plt.show()

def fof_scatter(halo,**kwargs):
    if 'snap_num' not in kwargs:
        if 'redshift' in kwargs:
            target_redshift=kwargs['redshift']
            snap_num,snap_redshift=get_snap_num(halo,target_redshift)
        else:
            sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    else:
        snap_num=kwargs['snap_num']
        snap_redshift=get_redshift(halo,snap_num)

    sf_positions=np.load(f'halos/{halo}/{snap_num}/subfind/fof_positions.npy')*units.Mpc
    print('\nSubfind Position Data Imported')
    sf_masses=np.load(f'halos/{halo}/{snap_num}/subfind/fof_masses.npy')*10**10*units.M_sun
    print('Subfind Mass Data Imported')
    
    fig_fofscatter,ax_fofscatter=plt.subplots(1,3,figsize=(15,15),constrained_layout=True)

    plane_indexes=[[0,1],[0,2],[1,2]]
    scatters=[]

    for plane in plane_indexes:
        scatter_axis=plane_indexes.index(plane)
        
        ax_fofscatter[scatter_axis].set_facecolor('black')
        scatter=ax_fofscatter[scatter_axis].scatter(sf_positions[:,plane[0]],sf_positions[:,plane[1]],s=.5,c=np.log10(sf_masses.to_value(units.kg)),cmap=cm.afmhot,vmax=34)
        scatters.append(scatter)
        
        fof1_marker=ax_fofscatter[scatter_axis].scatter(sf_positions[0,plane[0]],sf_positions[0,plane[1]],c='r',marker='*',s=100)

    colourbar=fig_fofscatter.colorbar(scatters[-1],ax=ax_fofscatter,shrink=.25)
    colourbar.set_label('$log_{10}($FoF Group Mass$)$ ($kg$)',fontsize=10)

    
    ax_fofscatter[0].set_xlabel('$x$ ($Mpc$)')
    ax_fofscatter[0].set_ylabel('$y$ ($Mpc$)')


    ax_fofscatter[1].set_xlabel('$x$ ($Mpc$)')
    ax_fofscatter[1].set_ylabel('$z$ ($Mpc$)')

    display_redshift=np.round(get_redshift(halo,snap_num),3)
    display_halo=halo.replace('_',' ')

    ax_fofscatter[1].set_title(f'{display_halo}, $z = {display_redshift}$')


    ax_fofscatter[2].set_xlabel('$y$ ($Mpc$)')
    ax_fofscatter[2].set_ylabel('$z$ ($Mpc$)')
    
    for ax in ax_fofscatter:
            ax.set_box_aspect(1.0)
            ax.xaxis.label.set_size(12)
            ax.yaxis.label.set_size(12)
            ax.tick_params(labelsize=10)

    plt.show()


def proj_mass_density(halo,bin_num,**kwargs):
    if 'snap_num' not in kwargs:
        if 'redshift' in kwargs:
            target_redshift=kwargs['redshift']
            snap_num,snap_redshift=get_snap_num(halo,target_redshift)
        else:
            sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    else:
        snap_num=kwargs['snap_num']
        snap_redshift=get_redshift(halo,snap_num)

    display_redshift=np.round(snap_redshift,3)
    display_halo=halo.replace('_',' ')

    if 'type_only' in kwargs:
        req_types=[kwargs['type_only']]
    else:
        req_types=['gas','dm','stars']

    type_plot_info={'gas':{'cmap':'plasma','title':f'{display_halo} Gas, z =${display_redshift}$, {bin_num} bins'},'dm':{'cmap':'viridis','title':f'{display_halo} Dark Matter, z =${display_redshift}$, {bin_num} bins'},'stars':{'cmap':'magma','title':f'{display_halo} Stars, z =${display_redshift}$, {bin_num} bins'}}

    fig_masshist,ax_masshist=plt.subplots(len(req_types),3,figsize=(15,15),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    for matter_type in req_types:
        loaded_data={f'{plane}':np.load(f'halos/{halo}/{snap_num}/binned/{bin_num}px/total_mass/{matter_type}.npy')[planes[plane]['index']]for plane in planes}
        loaded_data['rel_pos']=read_raw_file(halo,matter_type,'rel_pos',snap_num=snap_num)
        
        vmin=np.log10(np.min([np.min(loaded_data[plane][loaded_data[plane]!=0.0]) for plane in planes])) #Find minimum (excluding zeros) and maximum projected density values across projection axes to normalise colour bars to
        vmax=np.log10(np.max([loaded_data[plane]for plane in planes]))

        extents=calc.get_extent(loaded_data['rel_pos'])
        plot_extents={plane:[extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc), extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)] for plane in planes}

        row_index=req_types.index(matter_type)

        imshows=[ax_masshist[row_index][planes[plane]['index']].imshow(np.log10(loaded_data[plane]),extent=plot_extents[plane],vmin=vmin,vmax=vmax,cmap=type_plot_info[matter_type]['cmap'],aspect='equal') for plane in planes]
        for plane in planes:
            ax_masshist[row_index][planes[plane]['index']].set_xlabel(planes[plane]['x_label'])
            ax_masshist[row_index][planes[plane]['index']].set_ylabel(planes[plane]['y_label'])
            if planes[plane]['index']==1:
                ax_masshist[row_index][planes[plane]['index']].set_title(type_plot_info[matter_type]['title'],fontsize=16,pad=20)
        
        colourbar=fig_masshist.colorbar(imshows[-1],ax=ax_masshist[row_index],shrink=.75)
        colourbar.set_label('$log_{10}($Projected Density$)$ ($g/cm^2$)')
        
        for ax in ax_masshist[row_index]:
            ax.set_aspect('equal', adjustable='box')
            
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(12)
            ax.yaxis.label.set_size(12)
            ax.tick_params(labelsize=10)

    plt.show()

def radial_mass_density(halo,**kwargs):
    if 'snap_num' not in kwargs:
        if 'redshift' in kwargs:
            target_redshift=kwargs['redshift']
            snap_num,snap_redshift=get_snap_num(halo,target_redshift)
        else:
            sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    else:
        snap_num=kwargs['snap_num']
        snap_redshift=get_redshift(halo,snap_num)

    display_redshift=np.round(snap_redshift,3)
    display_halo=halo.replace('_',' ')

    if 'type_only' in kwargs:
        req_types=[kwargs['type_only']]
    else:
        req_types=['gas','dm','stars']

    halo_r200=read_subfind_params(halo,snap_num=snap_num)['halo_r200'].to_value(units.kpc)

    type_plot_info={'gas':{'colour':'red','label':'Gas','label_pad':0.045},'dm':{'colour':'blue','label':'DM','label_pad':0.045},'stars':{'colour':'orange','label':'Stars','label_pad':0.055}}

    fig_radialdensity, ax_radialdensity = plt.subplots()

    for matter_type in req_types:
        loaded_data=np.load(f'halos/{halo}/{snap_num}/binned/radial_mass_density/{matter_type}.npy')
        
        bin_centres=loaded_data[0]
        densities=loaded_data[1]

        ax_radialdensity.plot(bin_centres,densities,c=type_plot_info[matter_type]['colour'])
        plt.text(bin_centres[0]-type_plot_info[matter_type]['label_pad'], densities[0], type_plot_info[matter_type]['label'], fontsize=12,c=type_plot_info[matter_type]['colour'])

    ax_radialdensity.axvline(halo_r200,c='k',ls='dashed')
    plt.text(halo_r200+1, 10**-27, '$R_{200_{crit}}=$'+str(np.round(halo_r200,1))+'$ kpc$', fontsize=12,c='k',rotation='vertical')

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Radius ($kpc$)')
    plt.ylabel('Spherical Radial Density ($g/cm^3$)')
    plt.xlim(left=0.05)
    plt.ylim(top=10**-22)
    plt.title(display_halo+', $z=$'+str(display_redshift))
    plt.show()

        
        

