from imports import *
from halo_readers import get_redshift, get_snap_num
from npy_data_readers import read_raw_file, read_subfind_params
import processing as calc

def dz_snapshot(halo,snap_num,**redshifts):   
    snap_nums,snap_redshifts=np.float64(np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/redshifts.npy'))

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
    plt.text(int(snap_num)+2, 32, f'Snapshot Number $= {snap_num}$', fontsize=12,c='k',rotation='vertical')

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

        if os.path.isdir(f'figures/{halo}/dz_snapshot/targets')!=True:
            os.makedirs(f'figures/{halo}/dz_snapshot/targets')
        
        plt.savefig(f'figures/{halo}/dz_snapshot/targets/{snap_num}.pdf',format="pdf",dpi=250,bbox_inches='tight')
    else:
        plt.savefig(f'figures/{halo}/dz_snapshot/{snap_num}.pdf',format="pdf",dpi=250,bbox_inches='tight')

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

    sf_positions=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/subfind/fof_positions.npy')*units.Mpc
    print('\nSubfind Position Data Imported')
    sf_masses=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/subfind/fof_masses.npy')*10**10*units.M_sun
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
    colourbar.set_label('$log_{10}($FoF Group Mass$)$ ($kg$)',fontsize=16)

    
    ax_fofscatter[0].set_xlabel('$x$ ($Mpc$)')
    ax_fofscatter[0].set_ylabel('$y$ ($Mpc$)')


    ax_fofscatter[1].set_xlabel('$x$ ($Mpc$)')
    ax_fofscatter[1].set_ylabel('$z$ ($Mpc$)')

    display_redshift=np.round(get_redshift(halo,snap_num),3)
    display_halo=halo.replace('_',' ')

    ax_fofscatter[1].set_title(f'{display_halo}, $z = {display_redshift}$',fontsize=18)


    ax_fofscatter[2].set_xlabel('$y$ ($Mpc$)')
    ax_fofscatter[2].set_ylabel('$z$ ($Mpc$)')
    
    for ax in ax_fofscatter:
            ax.set_box_aspect(1.0)
            ax.xaxis.label.set_size(16)
            ax.yaxis.label.set_size(16)
            ax.tick_params(labelsize=14)

    plt.savefig(f'figures/{halo}/fof_scatter/{snap_num}.png',format="png",dpi=250,bbox_inches='tight')

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

    type_plot_info={'gas':{'cmap':'plasma','title':f'{display_halo} Gas, z =${display_redshift}$, {bin_num}px'},'dm':{'cmap':'viridis','title':f'{display_halo} Dark Matter, z =${display_redshift}$, {bin_num}px'},'stars':{'cmap':'magma','title':f'{display_halo} Stars, z =${display_redshift}$, {bin_num}px'}}

    fig_masshist,ax_masshist=plt.subplots(len(req_types),3,figsize=(15,15),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    for matter_type in req_types:
        loaded_data={f'{plane}':np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/total_mass/{matter_type}.npy')[planes[plane]['index']]for plane in planes}
        loaded_data['rel_pos']=read_raw_file(halo,matter_type,'rel_pos',snap_num=snap_num)
        
        print(matter_type)
        print(np.min(loaded_data['xy']))
        print(np.max(loaded_data['xy']))

        vmin=np.log10(np.min([np.min(loaded_data[plane][loaded_data[plane]!=0.0]) for plane in planes])) #Find minimum (excluding zeros) and maximum projected density values across projection axes to normalise colour bars to
        vmax=np.log10(np.max([loaded_data[plane]for plane in planes]))

        extents=calc.get_extent(loaded_data['rel_pos'])
        plot_extents={plane:[extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc), extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)] for plane in planes}

        row_index=req_types.index(matter_type)

        imshows=[ax_masshist[row_index][planes[plane]['index']].imshow(np.log10(loaded_data[plane]),extent=plot_extents[plane],vmin=vmin,vmax=vmax,cmap=type_plot_info[matter_type]['cmap'],aspect='equal') for plane in planes]
        for plane in planes:
            ax_masshist[row_index][planes[plane]['index']].set_xlabel(planes[plane]['x_label'],fontsize=18,labelpad=2)
            ax_masshist[row_index][planes[plane]['index']].set_ylabel(planes[plane]['y_label'],fontsize=18,labelpad=2)
            if planes[plane]['index']==1:
                ax_masshist[row_index][planes[plane]['index']].set_title(type_plot_info[matter_type]['title'],fontsize=20,pad=20)
        
        colourbar=fig_masshist.colorbar(imshows[-1],ax=ax_masshist[row_index],shrink=.85)
        colourbar.set_label('$log_{10}($Projected Density$)$ ($g/cm^2$)',fontsize=18)
        colourbar.ax.tick_params(labelsize=16)
        
        for ax in ax_masshist[row_index]:
            ax.set_aspect('equal', adjustable='box')
            
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(18)
            ax.yaxis.label.set_size(18)

            ax.yaxis.set_label_coords(-0.1, 0.5)

            ax.tick_params(labelsize=16)

    plt.savefig(f'figures/{halo}/{bin_num}/proj_mass_density/{snap_num}.pdf',format="pdf",dpi=250,bbox_inches='tight')

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
        loaded_data=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/radial_mass_density/{matter_type}.npy')
        
        bin_centres=loaded_data[0]
        densities=loaded_data[1]

        ax_radialdensity.plot(bin_centres,densities,c=type_plot_info[matter_type]['colour'])
        plt.text(bin_centres[0]-type_plot_info[matter_type]['label_pad'], densities[0], type_plot_info[matter_type]['label'], fontsize=14,c=type_plot_info[matter_type]['colour'])

    ax_radialdensity.axvline(halo_r200,c='k',ls='dashed')
    plt.text(halo_r200+1, 10**-27, '$R_{200_{crit}}=$'+str(np.round(halo_r200,1))+'$ kpc$', fontsize=14,c='k',rotation='vertical')

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Radius ($kpc$)',fontsize=14)
    plt.ylabel('Spherical Radial Density ($g/cm^3$)',fontsize=14)
    plt.xlim(left=0.05)
    plt.ylim(top=10**-22)
    plt.title(display_halo+', $z=$'+str(display_redshift),fontsize=16)
    ax_radialdensity.xaxis.set_tick_params(labelsize=12)
    ax_radialdensity.yaxis.set_tick_params(labelsize=12)


    plt.savefig(f'figures/{halo}/radial_mass_density/{snap_num}.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def proj_gas_densities(halo,bin_num,**kwargs):
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

    if 'dens_only' in kwargs:
        req_dens=[kwargs['dens_only']]
    else:
        req_dens=['total_mass','hydrogen_mass','nH_col']

    dens_plot_info={'total_mass':{'cmap':'plasma','title':f'{display_halo} Total Gas, z =${display_redshift}$, {bin_num} bins','cbar_label':'$log_{10}($Projected Density$)$ ($g/cm^2$)'},'hydrogen_mass':{'cmap':'viridis','title':f'{display_halo} Hydrogen Mass Fraction, z =${display_redshift}$, {bin_num} bins','cbar_label':'$log_{10}($Projected Density$)$ ($g/cm^2$)'},'nH_col':{'cmap':'magma','title':f'{display_halo} Neutral Hydrogen Fraction, z =${display_redshift}$, {bin_num} bins','cbar_label':'$log_{10}($Projected Column Density$)$ ($n_{H_1^1}/cm^2$)'}}

    fig_gashist,ax_gashist=plt.subplots(len(req_dens),3,figsize=(15,15),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    for density in req_dens:
        loaded_data={f'{plane}':np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/{density}.npy')[planes[plane]['index']]for plane in planes}
        loaded_data['rel_pos']=read_raw_file(halo,'gas','rel_pos',snap_num=snap_num)
        
        vmin=np.log10(np.min([np.min(loaded_data[plane][loaded_data[plane]!=0.0]) for plane in planes])) #Find minimum (excluding zeros) and maximum projected density values across projection axes to normalise colour bars to
        vmax=np.log10(np.max([loaded_data[plane]for plane in planes]))

        extents=calc.get_extent(loaded_data['rel_pos'])
        plot_extents={plane:[extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc), extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)] for plane in planes}

        row_index=req_dens.index(density)

        imshows=[ax_gashist[row_index][planes[plane]['index']].imshow(np.log10(loaded_data[plane]),vmin=vmin,vmax=vmax,extent=plot_extents[plane],cmap=dens_plot_info[density]['cmap'],aspect='equal') for plane in planes]
        for plane in planes:
            ax_gashist[row_index][planes[plane]['index']].set_xlabel(planes[plane]['x_label'])
            ax_gashist[row_index][planes[plane]['index']].set_ylabel(planes[plane]['y_label'])
            if planes[plane]['index']==1:
                ax_gashist[row_index][planes[plane]['index']].set_title(dens_plot_info[density]['title'],fontsize=20,pad=20)
        
        colourbar=fig_gashist.colorbar(imshows[-1],ax=ax_gashist[row_index],shrink=.9)
        colourbar.set_label(dens_plot_info[density]['cbar_label'],fontsize=18)
        colourbar.ax.tick_params(labelsize=16)
        
        for ax in ax_gashist[row_index]:
            ax.set_aspect('equal', adjustable='box')
            
            ax.set_aspect('equal', adjustable='box')
            
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(18)
            ax.yaxis.label.set_size(18)

            ax.yaxis.set_label_coords(-0.1, 0.5)

            ax.tick_params(labelsize=16)

    plt.savefig(f'figures/{halo}/{bin_num}/proj_gas_densities/{snap_num}.pdf',format="pdf",dpi=250,bbox_inches='tight')
    
    plt.show()

def weighted_mean_gz(halo,bin_num,**kwargs): 
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

    fig_gzhist,ax_gzhist=plt.subplots(1,3,figsize=(15,15),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    loaded_data={f'{plane}':np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/mean_gz.npy')[planes[plane]['index']]for plane in planes}
    loaded_data['rel_pos']=read_raw_file(halo,'gas','rel_pos',snap_num=snap_num)
    
    all_planes_data=np.array([loaded_data['xy'],loaded_data['xz'],loaded_data['yz']])
    vmin=np.log10(np.nanmin(all_planes_data)) 
    vmax=np.log10(np.nanmax(all_planes_data))

    extents=calc.get_extent(loaded_data['rel_pos'])
    plot_extents={plane:[extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc), extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)] for plane in planes}

    imshows=[ax_gzhist[planes[plane]['index']].imshow(np.log10(loaded_data[plane]),extent=plot_extents[plane],vmin=vmin,vmax=vmax,cmap='plasma',aspect='equal') for plane in planes]

    for plane in planes:
            ax_gzhist[planes[plane]['index']].set_xlabel(planes[plane]['x_label'])
            ax_gzhist[planes[plane]['index']].set_ylabel(planes[plane]['y_label'])
            if planes[plane]['index']==1:
                ax_gzhist[planes[plane]['index']].set_title(f'{display_halo} Mass-Weighted Mean Metallicity, z =${display_redshift}$, {bin_num} bins',fontsize=20,pad=20)
        
    colourbar=fig_gzhist.colorbar(imshows[-1],ax=ax_gzhist,shrink=.25)
    colourbar.set_label('$log_{10}($Mass-Weighted Mean Solar-Relative Metallicity$)$ ($Z_\odot$)',fontsize=12)
    colourbar.ax.tick_params(labelsize=16)
        
    for ax in ax_gzhist:
        ax.set_aspect('equal', adjustable='box')
        
        ax.set_box_aspect(1.0)
    
        ax.xaxis.label.set_size(18)
        ax.yaxis.label.set_size(18)

        ax.yaxis.set_label_coords(-0.1, 0.5)

        ax.tick_params(labelsize=16)

    plt.savefig(f'figures/{halo}/{bin_num}/weighted_mean_gz/{snap_num}.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def nH_col_gz_scatter(halo,bin_num,plane,**kwargs):
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

    snap_num=str(snap_num)
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
        snap_num='0'+snap_num


    if 'scatter_only' in kwargs:
        scatter_only=kwargs['scatter_only']
    else:
        scatter_only=False

    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    if plane not in planes:
        sys.exit('Please provide a cartesian plane (\"plane=ab\")')

    halo_r200=read_subfind_params(halo,snap_num=snap_num)['halo_r200'].value
    
    nH_col=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/nH_col.npy')[planes[plane]['index']].flatten()
    mean_gz=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/mean_gz.npy')[planes[plane]['index']].flatten()
    bin_radii=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/bin_radii.npy')[planes[plane]['index']].flatten()/halo_r200

    masked_data={'DLA':{'param':'mean_gz'},'subDLA':{'param':'mean_gz'},'LymanLimit':{'param':'mean_gz'},'lo_z':{'param':'nH_col'}}
    
    for mask in masked_data:
        loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked/{mask}.npz')
        
        masked_data[mask]['data'] = np.ma.masked_array(loaded[f'{masked_data[mask]["param"]}_data'], mask=loaded[f'{masked_data[mask]["param"]}_mask'])

    lo_z_nH_col=masked_data['lo_z']['data'][planes[plane]['index']].compressed()
    DLA_mean_gz=masked_data['DLA']['data'][planes[plane]['index']].compressed()
    subDLA_mean_gz=masked_data['subDLA']['data'][planes[plane]['index']].compressed()
    LymanLimit_mean_gz=masked_data['LymanLimit']['data'][planes[plane]['index']].compressed()

    fig_scatter, ax_scatter=plt.subplots(figsize=(15,12))

    scatter=ax_scatter.scatter(nH_col,mean_gz,marker='x',c=bin_radii,cmap='plasma_r',s=2,alpha=.2,zorder=10)

    plt.xlabel(f'Pixel Projected ${plane}$ Planar Number Density ($H_1^1/cm^2$)',fontsize=18)
    plt.ylabel('Solar-Relative Pixel-Mass-Weighted Mean Metallicity ($Z_\odot$)',fontsize=18)
    ax_scatter.tick_params(labelsize=16)

    scatter_colourbar=fig_scatter.colorbar(scatter, ax=ax_scatter,location='left')
    scatter_colourbar.solids.set_alpha(1)
    scatter_colourbar.ax.tick_params(labelsize=16)
    scatter_colourbar.set_label('$R_{200_{crit}}$-Normalised Radial Distance From Centre Of FoF Group ($R_{200_{crit}}$)',fontsize=20)

    xlims=[np.float64(10**9),np.float64(10**23)]
    ylims=[np.float64(10**-6),np.float64(10**1)]
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(xlims)
    plt.ylim(ylims)

    ax_scatter.axvline(np.float64(10**20.3),c='r',ls='dashed')
    ax_scatter.fill_betweenx(np.array([10**-8,10**4]),np.float64(10**20.3),np.float64(10**23),color='r',alpha=.2)
    ax_scatter.text(np.float64(0.2*10**22),np.float64(10**-2.5),'DLA',c='r',rotation=45,fontsize=14)

    ax_scatter.axvline(10**19,c='b',ls='dashed')
    ax_scatter.fill_betweenx(np.array([10**-8,10**4]),np.float64(10**19),np.float64(10**20.3),color='b',alpha=.2)
    ax_scatter.text(np.float64(1.3*10**19),np.float64(10**-2.5),'Sub-DLA',c='b',rotation=45,fontsize=14)
    
    ax_scatter.axvline(10**17.2,c='g',ls='dashed')
    ax_scatter.fill_betweenx(np.array([10**-8,10**4]),np.float64(10**17.2),np.float64(10**19),color='g',alpha=.2)
    ax_scatter.text(np.float64(2.1*10**17),np.float64(10**-2.5),'Lyman Limit',c='g',rotation=45,fontsize=14)

    ax_scatter.axhline(10**-3,color='blueviolet',ls='dashed',lw=2)
    ax_scatter.text(np.float64(0.8*10**14),np.float64(1.2*10**-3),'Low Metallicity Threshold, $Z \leq 10^{-3}Z_{\odot}$',fontsize=18,color='blueviolet')

    ax_scatter.axhline(1,color='indigo',ls='dashed',lw=2)
    ax_scatter.text(np.float64(0.8*10**16),np.float64(1.2),'Solar Metallicity, $Z \leq Z_{\odot}$',fontsize=18,color='indigo')
    
    if scatter_only==True:
        if os.path.isdir(f'figures/{halo}/{bin_num}/nH_col_gz_scatter/scatter_only')!=True:
            os.makedirs(f'figures/{halo}/{bin_num}/nH_col_gz_scatter/scatter_only')
        plt.savefig(f'figures/{halo}/{bin_num}/nH_col_gz_scatter/scatter_only/{snap_num}.png',format="png",dpi=250,bbox_inches='tight')

    if scatter_only==False:
        left,bottom,width,height=ax_scatter.get_position().bounds

        axis_bins=100

        ax_projdenshist=fig_scatter.add_axes([left,bottom+height,width,height/4],sharex=ax_scatter)
        
        plt.hist(nH_col,bins=np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),axis_bins),zorder=10,color='white',edgecolor='black',log=True)
        plt.hist(lo_z_nH_col,bins=np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),axis_bins),zorder=11,color='blueviolet',edgecolor='purple',alpha=.5,log=True)
        
        plt.yscale('log')
        plt.ylim([0,np.float64(10**5)])
        plt.title(f'{display_halo} {plane} Projection, z =${display_redshift}$, {bin_num} bins',fontsize=20,pad=20)

        ax_projdenshist.tick_params(labelbottom=False,labelleft=False,labelright=True,labelsize=16)
        ax_projdenshist.yaxis.tick_right()
        ax_projdenshist.set_ylabel('Number Density',fontsize=18)

        ax_projdenshist.axvline(np.float64(10**20.3),c='r',ls='dashed',zorder=11)
        ax_projdenshist.fill_betweenx([np.float64(0),np.float64(10**6)],np.float64(10**20.3),np.float64(10**23),color='r',alpha=.2)

        ax_projdenshist.axvline(10**19,c='b',ls='dashed',zorder=11)
        ax_projdenshist.fill_betweenx([np.float64(0),np.float64(10**6)],np.float64(10**19),np.float64(10**20.3),color='b',alpha=.2)

        ax_projdenshist.axvline(10**17.2,c='g',ls='dashed',zorder=11)
        ax_projdenshist.fill_betweenx([np.float64(0),np.float64(10**6)],np.float64(10**17.2),np.float64(10**19),color='g',alpha=.2,zorder=1)


        ax_gzhist=fig_scatter.add_axes([left+width,bottom,width/4,height],sharey=ax_scatter)
        plt.hist(mean_gz,bins=np.logspace(np.log10(ylims[0]),np.log10(ylims[1]),axis_bins),orientation='horizontal',color='white',edgecolor='black',log=True)
        
        plt.hist(subDLA_mean_gz,bins=np.logspace(np.log10(ylims[0]),np.log10(ylims[1]),axis_bins),orientation='horizontal',color='blue',edgecolor='cyan',alpha=.6,log=True,lw=2)
        plt.hist(LymanLimit_mean_gz,bins=np.logspace(np.log10(ylims[0]),np.log10(ylims[1]),axis_bins),orientation='horizontal',color='green',edgecolor='springgreen',alpha=.5,log=True,lw=2)
        plt.hist(DLA_mean_gz,bins=np.logspace(np.log10(ylims[0]),np.log10(ylims[1]),axis_bins),orientation='horizontal',color='red',edgecolor='orangered',alpha=.4,log=True,lw=2)

        ax_gzhist.axhline(10**-3,color='blueviolet',ls='dashed',lw=2)
        ax_gzhist.axhline(1,color='indigo',ls='dashed',lw=2)
        ax_gzhist.set_xlim([0,np.float64(3*10**4)])
        ax_gzhist.set_xscale('log')

        ax_gzhist.tick_params(labelleft=False,labelbottom=False,labeltop=True,labelsize=16)
        ax_gzhist.xaxis.tick_top()
        ax_gzhist.set_xlabel('Number Density',fontsize=18)
        gzhist_labels=ax_gzhist.get_xticklabels()
        gzhist_labels[0].set_visible(False)

        plt.savefig(f'figures/{halo}/{bin_num}/nH_col_gz_scatter/{snap_num}/{plane}.png',format="png",dpi=250,bbox_inches='tight')

    plt.show()


def threshold_mass_fracs(halo,bin_num,plane,**kwargs):
    if os.path.isdir(f'figures/{halo}/{bin_num}/masked_vs_redshifts')!=True:
        os.makedirs(f'figures/{halo}/{bin_num}/masked_vs_redshifts/mass_fracs')

    display_halo=halo.replace('_',' ')

    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}
    if plane not in planes:
        sys.exit('Please provide a cartesian plane (\"plane=ab\")')
    
    plane_index=planes[plane]['index']

    redshifts=np.load(f'halos/{halo}/mass_fracs/{bin_num}px/redshifts.npy')
    mass_frac_DLAs=np.load(f'halos/{halo}/mass_fracs/{bin_num}px/DLA.npy')[plane_index]  
    mass_frac_subDLAs=np.load(f'halos/{halo}/mass_fracs/{bin_num}px/subDLA.npy')[plane_index]
    mass_frac_LymanLimits=np.load(f'halos/{halo}/mass_fracs/{bin_num}px/LymanLimit.npy')[plane_index]
    mass_frac_lo_z_DLAs=np.load(f'halos/{halo}/mass_fracs/{bin_num}px/lo_z_DLA.npy')[plane_index]  
    mass_frac_lo_z_subDLAs=np.load(f'halos/{halo}/mass_fracs/{bin_num}px/lo_z_subDLA.npy')[plane_index]
    mass_frac_lo_z_LymanLimits=np.load(f'halos/{halo}/mass_fracs/{bin_num}px/lo_z_LymanLimit.npy')[plane_index]

    fig_time,ax_time=plt.subplots()

    ax_time.scatter(redshifts,np.array(mass_frac_DLAs),c='r',marker='d')
    ax_time.scatter(redshifts,np.array(mass_frac_subDLAs),c='b',marker='d')
    ax_time.scatter(redshifts,np.array(mass_frac_LymanLimits),c='g',marker='d')

    ax_time.scatter(redshifts,np.array(mass_frac_lo_z_DLAs),c='w',edgecolor='r',marker='d')
    ax_time.vlines(redshifts,np.array(mass_frac_DLAs),np.array(mass_frac_lo_z_DLAs),colors='r',ls='dashed',zorder=0,alpha=.5)

    ax_time.scatter(redshifts,np.array(mass_frac_lo_z_subDLAs),c='w',edgecolor='b',marker='d')
    ax_time.vlines(redshifts,np.array(mass_frac_subDLAs),np.array(mass_frac_lo_z_subDLAs),colors='b',ls='dashed',zorder=0,alpha=.5)

    ax_time.scatter(redshifts,np.array(mass_frac_lo_z_LymanLimits),c='w',edgecolor='g',marker='d')
    ax_time.vlines(redshifts,np.array(mass_frac_LymanLimits),np.array(mass_frac_lo_z_LymanLimits),colors='g',ls='dashed',zorder=0,alpha=.5)
    
    ax_time.invert_xaxis()

    ax_time.xaxis.label.set_size(18)
    ax_time.yaxis.label.set_size(18)
    ax_time.tick_params(labelsize=18)
    plt.title(f'{display_halo}, {plane}-plane, {bin_num} bins',fontsize=20)
    plt.yscale('log')
    plt.xlabel('Redshift (z)',fontsize=18)
    plt.ylabel(r'Mass Fraction ($\frac{m_{px_{threshold}}}{M_{halo}}$)',fontsize=18)

    plt.savefig(f'figures/{halo}/{bin_num}/masked_vs_redshifts/mass_fracs/{plane}.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def contour_gas_hists(halo,bin_num,**kwargs):
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
    
    snap_num=str(snap_num)
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
        snap_num='0'+snap_num

    dens_plot_info={'cmap':'plasma','title':f'{display_halo} Total Gas'}

    fig_gashist,ax_gashist=plt.subplots(1,3,figsize=(15,15),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}
    
    loaded_data={f'{plane}':np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/total_mass.npy')[planes[plane]['index']]for plane in planes}
    loaded_data['rel_pos']=read_raw_file(halo,'gas','rel_pos',snap_num=snap_num)

    masks={'DLA':{'colour':'r','cmap':ListedColormap(np.array([[1,0,0,.6],[1,0,0,0]]))},'subDLA':{'colour':'b','cmap':ListedColormap(np.array([[0,0,1,.3],[0,0,1,0]]))},'LymanLimit':{'colour':'g','cmap':ListedColormap(np.array([[0,1,0,.3],[0,1,0,0]]))}}
    
    for mask in masks:
        loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked/{mask}.npz') 
        masks[mask]['data'] = np.ma.masked_array(loaded[f'nH_col_data'], mask=loaded[f'nH_col_mask'])

        lo_z_loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked/lo_z_{mask}.npz')
        masks[mask]['lo_z_data'] = np.ma.masked_array(lo_z_loaded[f'nH_col_data'], mask=lo_z_loaded[f'nH_col_mask'])
        
    vmin=np.log10(np.min([np.min(loaded_data[plane][loaded_data[plane]!=0.0]) for plane in planes])) #Find minimum (excluding zeros) and maximum projected density values across projection axes to normalise colour bars to
    vmax=np.log10(np.max([loaded_data[plane]for plane in planes]))

    extents=calc.get_extent(loaded_data['rel_pos'])
    plot_extents={plane:[extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc), extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)] for plane in planes}

    imshows=[ax_gashist[planes[plane]['index']].imshow(np.log10(loaded_data[plane]),extent=plot_extents[plane],vmin=vmin,vmax=vmax,cmap=dens_plot_info['cmap'],aspect='equal') for plane in planes]

    LymanLimit_contours=[ax_gashist[planes[plane]['index']].contour(masks['LymanLimit']['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=plot_extents[plane],vmin=vmin,vmax=vmax,colors=masks['LymanLimit']['colour']) for plane in planes if np.any(~masks['LymanLimit']['data'][planes[plane]['index']].mask)]
    LymanLimit_fills=[ax_gashist[planes[plane]['index']].imshow(masks['LymanLimit']['data'][planes[plane]['index']].mask.astype(float),extent=plot_extents[plane],cmap=masks['LymanLimit']['cmap']) for plane in planes if np.any(~masks['LymanLimit']['data'][planes[plane]['index']].mask)]
    lo_z_LymanLimit_contours=[ax_gashist[planes[plane]['index']].contourf(masks['LymanLimit']['lo_z_data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0,0.5],extent=plot_extents[plane],vmin=vmin,vmax=vmax,colors='none',hatches=['xxxx']) for plane in planes if np.any(~masks['LymanLimit']['lo_z_data'][planes[plane]['index']].mask)]
    for contour in lo_z_LymanLimit_contours:
        contour.set_edgecolor('blueviolet')
    
    subDLA_contours=[ax_gashist[planes[plane]['index']].contour(masks['subDLA']['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=plot_extents[plane],vmin=vmin,vmax=vmax,colors=masks['subDLA']['colour']) for plane in planes if np.any(~masks['subDLA']['data'][planes[plane]['index']].mask)]
    sub_DLA_fills=[ax_gashist[planes[plane]['index']].imshow(masks['subDLA']['data'][planes[plane]['index']].mask.astype(float),extent=plot_extents[plane],cmap=masks['subDLA']['cmap']) for plane in planes if np.any(~masks['subDLA']['data'][planes[plane]['index']].mask)]
    lo_z_subDLA_contours=[ax_gashist[planes[plane]['index']].contourf(masks['subDLA']['lo_z_data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0,0.5],extent=plot_extents[plane],vmin=vmin,vmax=vmax,colors='none',hatches=['xxxx']) for plane in planes if np.any(~masks['subDLA']['lo_z_data'][planes[plane]['index']].mask)]
    for contour in lo_z_subDLA_contours:
        contour.set_edgecolor('blueviolet')
    
    DLA_contours=[ax_gashist[planes[plane]['index']].contour(masks['DLA']['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=plot_extents[plane],vmin=vmin,vmax=vmax,colors=masks['DLA']['colour']) for plane in planes if np.any(~masks['DLA']['data'][planes[plane]['index']].mask)]
    DLA_fills=[ax_gashist[planes[plane]['index']].imshow(masks['DLA']['data'][planes[plane]['index']].mask.astype(float),extent=plot_extents[plane],cmap=masks['DLA']['cmap']) for plane in planes if np.any(~masks['DLA']['data'][planes[plane]['index']].mask)]  
    lo_z_DLA_contours=[ax_gashist[planes[plane]['index']].contourf(masks['DLA']['lo_z_data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0,0.5],extent=plot_extents[plane],vmin=vmin,vmax=vmax,colors='none',hatches=['xxxx']) for plane in planes if np.any(~masks['DLA']['lo_z_data'][planes[plane]['index']].mask)]
    for contour in lo_z_DLA_contours:
        contour.set_edgecolor('blueviolet')

    for plane in planes:
        ax_gashist[planes[plane]['index']].set_xlabel(planes[plane]['x_label'])
        ax_gashist[planes[plane]['index']].set_ylabel(planes[plane]['y_label'])
        if planes[plane]['index']==1:
            ax_gashist[planes[plane]['index']].set_title(dens_plot_info['title'],fontsize=20,pad=20)
    
    colourbar=fig_gashist.colorbar(imshows[-1],ax=ax_gashist,shrink=.25)
    colourbar.set_label('$log_{10}($Projected Density$)$ ($g/cm^2$)',fontsize=18)
    colourbar.ax.tick_params(labelsize=16)
    
    for ax in ax_gashist:
        ax.set_aspect('equal', adjustable='box')
        
        ax.set_aspect('equal', adjustable='box')
        
        ax.set_box_aspect(1.0)
    
        ax.xaxis.label.set_size(18)
        ax.yaxis.label.set_size(18)

        ax.yaxis.set_label_coords(-0.1, 0.5)

        ax.tick_params(labelsize=16)

    if os.path.isdir(f'figures/{halo}/{bin_num}/contour_gas_densities')!=True:
        os.makedirs(f'figures/{halo}/{bin_num}/contour_gas_densities',exist_ok=True)

    plt.savefig(f'figures/{halo}/{bin_num}/contour_gas_densities/{snap_num}.pdf',format="pdf",dpi=250,bbox_inches='tight')
    
    plt.show()


def stellar_masses_redshift(halos):
    if os.path.isdir(f'figures/all_halos/stellar_masses')!=True:
        os.makedirs(f'figures/all_halos/stellar_masses')

    fig_masses,ax_masses=plt.subplots()

    colours=['r','b','g']

    dla_ranges=[(0,5),(0,0),(0,0)]

    for halo in halos:
        display_halo=halo.replace('_',' ')

        redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
        masses=np.load(f'halos/{halo}/stellar_masses/stellar_masses.npy')
          
        ax_masses.plot(redshifts,masses,c=colours[halos.index(halo)])

        ax_masses.scatter(redshifts[dla_ranges[halos.index(halo)][0]:dla_ranges[halos.index(halo)][1]],masses[dla_ranges[halos.index(halo)][0]:dla_ranges[halos.index(halo)][1]],c=colours[halos.index(halo)],marker='d')

        if halo=='halo8':
            shift=0.4
        else:
            shift=0.6
        plt.text(redshifts[0]+shift,masses[0],display_halo,fontsize=14,c=colours[halos.index(halo)])
    

    ax_masses.invert_xaxis()

    ax_masses.xaxis.label.set_size(18)
    ax_masses.yaxis.label.set_size(18)
    ax_masses.tick_params(labelsize=18)
    plt.title(f'Total Halo Stellar Mass',fontsize=20)
    plt.xlabel('Redshift (z)',fontsize=18)
    plt.ylabel(r'Total Stellar Mass ($M_{\odot}$)',fontsize=18)
    plt.xlim([4.7,0.9])
    plt.yscale('log')

    plt.savefig(f'figures/all_halos/stellar_masses/stellar_masses.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def metallicity_radius(halo,**kwargs):
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

    radii=read_raw_file(halo,'gas','radii',snap_num=snap_num)
    metallicities=read_raw_file(halo,'gas','gz',snap_num=snap_num)
    masses=read_raw_file(halo,'gas','mass',snap_num=snap_num)

    fig_gzrad,ax_gzrad=plt.subplots()

    scatter=ax_gzrad.scatter(radii,metallicities,marker='x',s=2,c=np.log(masses.value),cmap='magma_r')

    colourbar=fig_gzrad.colorbar(scatter, ax=ax_gzrad,location='right',alpha=1)
    colourbar.set_label('log(Gas Particle Mass)')

    plt.yscale('log')
    plt.xlabel('Radial Distance From Center ($Mpc$)')
    plt.ylabel('Gas Particle Metallicity ($Z_\odot$)')

    plt.show()

def metallicities_redshift(halos):
    if os.path.isdir(f'figures/all_halos/metallicities_redshift')!=True:
        os.makedirs(f'figures/all_halos/metallicities_redshift')

    fig_zs,ax_zs=plt.subplots()

    colours=['r','b','g']

    for halo in halos:
        display_halo=halo.replace('_',' ')

        redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
        gas_zs=np.load(f'halos/{halo}/metallicity_behaviour/gas.npy')
        star_zs=np.load(f'halos/{halo}/metallicity_behaviour/stars.npy')

        print(gas_zs)
          
        ax_zs.plot(redshifts,gas_zs,c=colours[halos.index(halo)],linestyle='dashed')
        ax_zs.plot(redshifts,star_zs,c=colours[halos.index(halo)],linestyle='dotted')
        if halo=='halo8':
            shift=0.4
        else:
            shift=0.6
        plt.text(redshifts[0]+shift,gas_zs[0],f'{display_halo} Gas',fontsize=14,c=colours[halos.index(halo)])
        plt.text(redshifts[0]+shift,star_zs[0],f'{display_halo} Stars',fontsize=14,c=colours[halos.index(halo)])

    ax_zs.invert_xaxis()

    ax_zs.xaxis.label.set_size(18)
    ax_zs.yaxis.label.set_size(18)
    ax_zs.tick_params(labelsize=18)
    plt.title(f'Stellar and Gas Solar-Relative Metallicity',fontsize=20)
    plt.xlabel('Redshift (z)',fontsize=18)
    plt.ylabel(r'Metallicity ($Z_\odot$)',fontsize=18)
    plt.xlim([4.7,0.9])

    plt.savefig(f'figures/all_halos/metallicities_redshift/metallicities.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()


def metallicities_stellar_masses(halo):
    if os.path.isdir(f'figures/{halo}/metallicities_stellar_masses')!=True:
        os.makedirs(f'figures/{halo}/metallicities_stellar_masses')

    fig_z_mass,ax_z_mass=plt.subplots()

    display_halo=halo.replace('_',' ')

    redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
    masses=np.load(f'halos/{halo}/stellar_masses/stellar_masses.npy')
    gas_zs=np.load(f'halos/{halo}/metallicity_behaviour/gas.npy')
    star_zs=np.load(f'halos/{halo}/metallicity_behaviour/stars.npy')

    vmin=min(redshifts)
    vmax=max(redshifts)

    gas_scatter=ax_z_mass.scatter(masses,gas_zs,c=redshifts,marker='o',cmap='Reds',vmin=vmin,vmax=vmax,zorder=1,edgecolors='black')
    stars_scatter=ax_z_mass.scatter(masses,star_zs,c=redshifts,marker='d',cmap='Reds',vmin=vmin,vmax=vmax,zorder=1,edgecolors='black')

    scatter_colourbar=fig_z_mass.colorbar(gas_scatter,ax=ax_z_mass)    

    ax_z_mass.xaxis.label.set_size(18)
    ax_z_mass.yaxis.label.set_size(18)
    ax_z_mass.tick_params(labelsize=18)

    ax_z_mass.set_xlabel('Stellar Mass ($M_\odot$)')
    ax_z_mass.set_ylabel('Metallicity ($Z_\odot$)')
    

    scatter_colourbar.set_label('Snapshot Redshift')

    plt.title(f'{display_halo} Metallicities and Stellar Masses')

    plt.savefig(f'figures/{halo}/metallicities_stellar_masses/metallicity_mass.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def threshold_area_fracs(halo,bin_num,plane,**kwargs):
    if os.path.isdir(f'figures/{halo}/{bin_num}/masked_vs_redshifts/area_fracs')!=True:
        os.makedirs(f'figures/{halo}/{bin_num}/masked_vs_redshifts/area_fracs')

    display_halo=halo.replace('_',' ')

    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}
    if plane not in planes:
        sys.exit('Please provide a cartesian plane (\"plane=ab\")')
    
    plane_index=planes[plane]['index']

    redshifts=np.load(f'halos/{halo}/area_fracs/{bin_num}px/redshifts.npy')
    area_frac_DLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/DLA.npy')[plane_index]  
    area_frac_subDLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/subDLA.npy')[plane_index]
    area_frac_LymanLimits=np.load(f'halos/{halo}/area_fracs/{bin_num}px/LymanLimit.npy')[plane_index]
    area_frac_lo_z_DLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/lo_z_DLA.npy')[plane_index]  
    area_frac_lo_z_subDLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/lo_z_subDLA.npy')[plane_index]
    area_frac_lo_z_LymanLimits=np.load(f'halos/{halo}/area_fracs/{bin_num}px/lo_z_LymanLimit.npy')[plane_index]

    fig_area_time,ax_area_time=plt.subplots()

    ax_area_time.scatter(redshifts,np.array(area_frac_DLAs),c='r',marker='d')
    ax_area_time.scatter(redshifts,np.array(area_frac_subDLAs),c='b',marker='d')
    ax_area_time.scatter(redshifts,np.array(area_frac_LymanLimits),c='g',marker='d')

    ax_area_time.scatter(redshifts,np.array(area_frac_lo_z_DLAs),c='w',edgecolor='r',marker='d')
    ax_area_time.vlines(redshifts,np.array(area_frac_DLAs),np.array(area_frac_lo_z_DLAs),colors='r',ls='dashed',zorder=0,alpha=.5)

    ax_area_time.scatter(redshifts,np.array(area_frac_lo_z_subDLAs),c='w',edgecolor='b',marker='d')
    ax_area_time.vlines(redshifts,np.array(area_frac_subDLAs),np.array(area_frac_lo_z_subDLAs),colors='b',ls='dashed',zorder=0,alpha=.5)

    ax_area_time.scatter(redshifts,np.array(area_frac_lo_z_LymanLimits),c='w',edgecolor='g',marker='d')
    ax_area_time.vlines(redshifts,np.array(area_frac_LymanLimits),np.array(area_frac_lo_z_LymanLimits),colors='g',ls='dashed',zorder=0,alpha=.5)
    
    ax_area_time.invert_xaxis()

    ax_area_time.xaxis.label.set_size(18)
    ax_area_time.yaxis.label.set_size(18)
    ax_area_time.tick_params(labelsize=18)
    plt.title(f'{display_halo}, {plane}-plane, {bin_num} bins',fontsize=20)
    plt.yscale('log')
    plt.xlabel('Redshift (z)',fontsize=18)
    plt.ylabel(r'Surface Area Fraction ($\frac{N_{px_{threshold}}}{N_{total}}$)',fontsize=18)

    plt.savefig(f'figures/{halo}/{bin_num}/masked_vs_redshifts/area_fracs/{plane}.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def m200_redshift(halos):
    if os.path.isdir(f'figures/all_halos/m200_redshift')!=True:
        os.makedirs(f'figures/all_halos/m200_redshift')

    fig_zs,ax_zs=plt.subplots()

    colours=['r','b','g']

    for halo in halos:
        display_halo=halo.replace('_',' ')

        redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
        m200s=np.load(f'halos/{halo}/m200s/m200.npy')
          
        ax_zs.plot(redshifts,m200s,c=colours[halos.index(halo)])
        if halo=='halo8':
            shift=0.4
        else:
            shift=0.6
        plt.text(redshifts[0]+shift,m200s[0],f'{display_halo}',fontsize=14,c=colours[halos.index(halo)])


    ax_zs.invert_xaxis()

    ax_zs.xaxis.label.set_size(18)
    ax_zs.yaxis.label.set_size(18)
    ax_zs.tick_params(labelsize=18)
    plt.title(r'$M_{200_{crit}}$',fontsize=20)
    plt.xlabel('Redshift (z)',fontsize=18)
    plt.ylabel(r'$M_{200_{crit}} (M_\odot)$',fontsize=18)
    plt.xlim([4.7,0.9])

    plt.savefig(f'figures/all_halos/m200_redshift/m200s.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()


def rho_gz_scatter(halo,bin_num,plane,**kwargs):
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

    snap_num=str(snap_num)
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
        snap_num='0'+snap_num


    if 'scatter_only' in kwargs:
        scatter_only=kwargs['scatter_only']
    else:
        scatter_only=False

    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    if plane not in planes:
        sys.exit('Please provide a cartesian plane (\"plane=ab\")')

    halo_r200=read_subfind_params(halo,snap_num=snap_num)['halo_r200'].value
    
    rho=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/total_mass.npy')[planes[plane]['index']].flatten()
    mean_gz=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/mean_gz.npy')[planes[plane]['index']].flatten()
    bin_radii=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/bin_radii.npy')[planes[plane]['index']].flatten()/halo_r200

    fig_scatter, ax_scatter=plt.subplots(figsize=(15,12))

    scatter=ax_scatter.scatter(rho,mean_gz,marker='x',c=bin_radii,cmap='plasma_r',s=2,alpha=.2,zorder=10)

    plt.xlabel(f' Projected ${plane}$ Gas Column Density ($g/cm^2$)',fontsize=18)
    plt.ylabel('Solar-Relative Pixel-Mass-Weighted Mean Metallicity ($Z_\odot$)',fontsize=18)
    ax_scatter.tick_params(labelsize=16)

    scatter_colourbar=fig_scatter.colorbar(scatter, ax=ax_scatter,location='left')
    scatter_colourbar.solids.set_alpha(1)
    scatter_colourbar.ax.tick_params(labelsize=16)
    scatter_colourbar.set_label('$R_{200_{crit}}$-Normalised Radial Distance From Centre Of FoF Group ($R_{200_{crit}}$)',fontsize=20)

    xlims=[np.float64(10**-9),np.float64(10**-1)]
    ylims=[np.float64(10**-6),np.float64(10**1)]
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(xlims)
    plt.ylim(ylims)

    ax_scatter.axhline(10**-3,color='blueviolet',ls='dashed',lw=2)
    ax_scatter.text(np.float64(0.8*10**-6),np.float64(1.2*10**-3),'Low Metallicity Threshold, $Z \leq 10^{-3}Z_{\odot}$',fontsize=18,color='blueviolet')

    ax_scatter.axhline(1,color='indigo',ls='dashed',lw=2)
    ax_scatter.text(np.float64(0.8*10**-6),np.float64(1.2),'Solar Metallicity, $Z \leq Z_{\odot}$',fontsize=18,color='indigo')
    
    
    left,bottom,width,height=ax_scatter.get_position().bounds

    axis_bins=100

    ax_projdenshist=fig_scatter.add_axes([left,bottom+height,width,height/4],sharex=ax_scatter)
    
    plt.hist(rho,bins=np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),axis_bins),zorder=10,color='white',edgecolor='black',log=True)
    
    plt.yscale('log')
    plt.ylim([0,np.float64(10**5)])
    plt.title(f'{display_halo} {plane} Projection, z =${display_redshift}$, {bin_num} bins',fontsize=20,pad=20)

    ax_projdenshist.tick_params(labelbottom=False,labelleft=False,labelright=True,labelsize=16)
    ax_projdenshist.yaxis.tick_right()
    ax_projdenshist.set_ylabel('Number Density',fontsize=18)


    ax_gzhist=fig_scatter.add_axes([left+width,bottom,width/4,height],sharey=ax_scatter)
    plt.hist(mean_gz,bins=np.logspace(np.log10(ylims[0]),np.log10(ylims[1]),axis_bins),orientation='horizontal',color='white',edgecolor='black',log=True)
    
    ax_gzhist.axhline(10**-3,color='blueviolet',ls='dashed',lw=2)
    ax_gzhist.axhline(1,color='indigo',ls='dashed',lw=2)
    ax_gzhist.set_xlim([0,np.float64(3*10**4)])
    ax_gzhist.set_xscale('log')

    ax_gzhist.tick_params(labelleft=False,labelbottom=False,labeltop=True,labelsize=16)
    ax_gzhist.xaxis.tick_top()
    ax_gzhist.set_xlabel('Number Density',fontsize=18)
    gzhist_labels=ax_gzhist.get_xticklabels()
    gzhist_labels[0].set_visible(False)

    if os.path.isdir(f'figures/{halo}/{bin_num}/rho_gz_scatter/{snap_num}') != True:
        os.makedirs(f'figures/{halo}/{bin_num}/rho_gz_scatter/{snap_num}',exist_ok=True)

    plt.savefig(f'figures/{halo}/{bin_num}/rho_gz_scatter/{snap_num}/{plane}.png',format="png",dpi=250,bbox_inches='tight')

    plt.show()

def weighting_hists(halo,bin_num,plane,**kwargs):
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

    snap_num=str(snap_num)
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
        snap_num='0'+snap_num

    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    weighted_gz=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/mean_gz.npy')[planes[plane]['index']].flatten()
    unweighted_gz=np.load(f'halos/{halo}/weighting_test/{snap_num}/unweighted_gz.npy')[planes[plane]['index']].flatten()

    fig_weighted,ax_weighted=plt.subplots()

    xlims=[np.float64(10**-6),np.float64(10**1)]

    plt.hist(unweighted_gz,color='white',edgecolor='black',log=True,bins=np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),100))
    plt.hist(weighted_gz,color='blue',edgecolor='blue',log=True,bins=np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),100),alpha=.4)
    
    plt.xlim(xlims)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Pixel Metallicity Value ($Z_\odot$)')
    plt.ylabel('Number Density')

    if os.path.isdir(f'figures/{halo}/weighting_testing') != True:
        os.makedirs(f'figures/{halo}/weighting_testing',exist_ok=True)

    plt.savefig(f'figures/{halo}/weighting_testing/weighting_test.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def sightlines_contours(halo,**kwargs):
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

    snap_num=str(snap_num)
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
        snap_num='0'+snap_num

    nh_masks=['DLA','subDLA','LymanLimit']

    fig_sightlineconts,ax_sightlineconts=plt.subplots(len(nh_masks),3,figsize=(15,15),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    dens_plot_info={'cmap':'plasma','title':f'{display_halo} Total Gas'}

    radii=['half','1','2','5']
    colours=['darkviolet','blue','dodgerblue','springgreen']

    loaded_data={f'{plane}':np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/512px/gas_only/total_mass.npy')[planes[plane]['index']]for plane in planes}

    vmin=np.log10(np.min([np.min(loaded_data[plane][loaded_data[plane]!=0.0]) for plane in planes])) #Find minimum (excluding zeros) and maximum projected density values across projection axes to normalise colour bars to
    vmax=np.log10(np.max([loaded_data[plane]for plane in planes]))

    obj_rel_pos=read_raw_file(halo,'gas','rel_pos',snap_num=snap_num)

    obj_extents=calc.get_extent(obj_rel_pos)
    raw_plot_extents={plane:np.array([obj_extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc),obj_extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc), obj_extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc),obj_extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)]) for plane in planes}

    for plane in planes:
        obj_span=[raw_plot_extents[plane][1]-raw_plot_extents[plane][0],raw_plot_extents[plane][3]-raw_plot_extents[plane][2]]
        planes[plane]['aspect_ratio']=obj_span[1]/obj_span[0]

    obj_plot_extents={plane:np.array([obj_extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc)*planes[plane]['aspect_ratio'],obj_extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc)*planes[plane]['aspect_ratio'], obj_extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc)/planes[plane]['aspect_ratio'],obj_extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)/planes[plane]['aspect_ratio']]) for plane in planes}

    for col_dens in nh_masks:
        for radius in radii:
            rad_masks={'DLA':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.6],[1,1,1,0]]))},'subDLA':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.3],[1,1,1,0]]))},'LymanLimit':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.3],[1,1,1,0]]))}}
            
            loaded = np.load(f'halos/{halo}/sightlines/{snap_num}/{radius}/{col_dens}.npz') 
            rad_masks[col_dens]['data'] = np.ma.masked_array(loaded[f'nH_col_data'], mask=loaded[f'nH_col_mask'])

            lo_z_loaded = np.load(f'halos/{halo}/sightlines/{snap_num}/{radius}/lo_z_{col_dens}.npz')
            rad_masks[col_dens]['lo_z_data'] = np.ma.masked_array(lo_z_loaded[f'nH_col_data'], mask=lo_z_loaded[f'nH_col_mask'])

            halo_pos=read_subfind_params(halo,snap_num=snap_num)['halo_pos']

            pos=np.load(f'halos/{halo}/sightlines/{snap_num}/{radius}/pos.npy')*units.Mpc
            rel_pos=pos-halo_pos

            extents=calc.get_extent(rel_pos)
            raw_extents={plane:np.array([extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc), extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)]) for plane in planes}
            
            raw_contours=[ax_sightlineconts[nh_masks.index(col_dens)][planes[plane]['index']].contour(rad_masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=raw_extents[plane],colors=colours[radii.index(radius)],zorder=3) for plane in planes if np.any(~rad_masks[col_dens]['data'][planes[plane]['index']].mask)]

            for contour in raw_contours:
                contour.set_visible(False)

            contour_centres={}

            for plane in planes:
                contour=raw_contours[planes[plane]['index']]
                shifted_extents=[]
                for path in contour.get_paths():
                    verts=path.vertices
                    centroid=np.mean(verts, axis=0)
                    contour_centres[plane]=centroid

            extents={plane:[raw_extents[plane][0]-contour_centres[plane][0],raw_extents[plane][1]-contour_centres[plane][0],raw_extents[plane][2]-contour_centres[plane][1],raw_extents[plane][3]-contour_centres[plane][1]] for plane in planes}

            contours=[ax_sightlineconts[nh_masks.index(col_dens)][planes[plane]['index']].contour(rad_masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=extents[plane],colors=colours[radii.index(radius)],zorder=3) for plane in planes if np.any(~rad_masks[col_dens]['data'][planes[plane]['index']].mask)]
            
        masks={'DLA':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.3],[1,1,1,0]]))},'subDLA':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.3],[1,1,1,0]]))},'LymanLimit':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.3],[1,1,1,0]]))}}

        loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/512px/gas_only/masked/{col_dens}.npz') 
        masks[col_dens]['data'] = np.ma.masked_array(loaded[f'nH_col_data'], mask=loaded[f'nH_col_mask'])

        lo_z_loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/512px/gas_only/masked/{col_dens}.npz')
        masks[col_dens]['lo_z_data'] = np.ma.masked_array(lo_z_loaded[f'nH_col_data'], mask=lo_z_loaded[f'nH_col_mask'])

        raw_obj_contours=[ax_sightlineconts[nh_masks.index(col_dens)][planes[plane]['index']].contour(masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=obj_plot_extents[plane],colors=masks[col_dens]['colour']) for plane in planes if np.any(~masks[col_dens]['data'][planes[plane]['index']].mask)]

        for contour in raw_obj_contours:
            contour.set_visible(False)

        obj_contour_centres={}

        for plane in planes:
            contour=raw_obj_contours[planes[plane]['index']]
            shifted_extents=[]
            for path in contour.get_paths():
                verts=path.vertices
                centroid=np.mean(verts, axis=0)
                obj_contour_centres[plane]=centroid

        shifted_obj_extents={plane:[obj_plot_extents[plane][0]-obj_contour_centres[plane][0],obj_plot_extents[plane][1]-obj_contour_centres[plane][0],obj_plot_extents[plane][2]-obj_contour_centres[plane][1],obj_plot_extents[plane][3]-obj_contour_centres[plane][1]] for plane in planes}

        adj_obj_contours=[ax_sightlineconts[nh_masks.index(col_dens)][planes[plane]['index']].contour(masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=shifted_obj_extents[plane],colors=masks[col_dens]['colour'],zorder=2) for plane in planes if np.any(~masks[col_dens]['data'][planes[plane]['index']].mask)]
        adj_obj_fills=[ax_sightlineconts[nh_masks.index(col_dens)][planes[plane]['index']].imshow(masks[col_dens]['data'][planes[plane]['index']].mask.astype(float),extent=shifted_obj_extents[plane],cmap=masks[col_dens]['cmap'],zorder=1) for plane in planes if np.any(~masks[col_dens]['data'][planes[plane]['index']].mask)]  
        adj_obj_imshows=[ax_sightlineconts[nh_masks.index(col_dens)][planes[plane]['index']].imshow(np.log10(loaded_data[plane]),extent=shifted_obj_extents[plane],vmin=vmin,vmax=vmax,cmap=dens_plot_info['cmap'],aspect='equal',zorder=0) for plane in planes]
       
        for plane in planes:
            ax_sightlineconts[nh_masks.index(col_dens)][planes[plane]['index']].set_xlabel(planes[plane]['x_label'])
            ax_sightlineconts[nh_masks.index(col_dens)][planes[plane]['index']].set_ylabel(planes[plane]['y_label'])
            if planes[plane]['index']==1:
                ax_sightlineconts[nh_masks.index(col_dens)][planes[plane]['index']].set_title(col_dens,fontsize=20,pad=20)
    
    radii=[0.5,1,2,5,10]

    row_dims=[1,1.5,2.5]
    row_index=0

    for row in ax_sightlineconts:
        for ax in row:
            ax.set_xlim(-row_dims[row_index],row_dims[row_index])
            ax.set_ylim(-row_dims[row_index],row_dims[row_index])
            ax.set_aspect('equal', adjustable='box')
            ax.set_aspect('equal', adjustable='box')
            ax.set_box_aspect(1.0)
            ax.xaxis.label.set_size(18)
            ax.yaxis.label.set_size(18)

            ax.yaxis.set_label_coords(-0.1, 0.5)

            ax.tick_params(labelsize=16)
        row_index+=1

    if os.path.isdir(f'figures/{halo}/sightlines/{snap_num}') != True:
        os.makedirs(f'figures/{halo}/sightlines/{snap_num}',exist_ok=True)

    plt.savefig(f'figures/{halo}/sightlines/{snap_num}/sightline_contours.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def sightlines_scatter(halo,plane,**kwargs):
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

    snap_num=str(snap_num)
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
        snap_num='0'+snap_num
        
    bin_num=512

    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    fig_sightlinescatter, ax_sightlinescatter=plt.subplots(figsize=(15,12))

    halo_r200=read_subfind_params(halo,snap_num=snap_num)['halo_r200'].value

    nH_col=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/nH_col.npy')[planes[plane]['index']].flatten()
    mean_gz=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/mean_gz.npy')[planes[plane]['index']].flatten()
    bin_radii=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/bin_radii.npy')[planes[plane]['index']].flatten()/halo_r200

    masked_data={'DLA':{'param':'mean_gz'},'subDLA':{'param':'mean_gz'},'LymanLimit':{'param':'mean_gz'}}

    for mask in masked_data:
        loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked/{mask}.npz')
        
        masked_data[mask]['data'] = np.ma.masked_array(loaded[f'{masked_data[mask]["param"]}_data'], mask=loaded[f'{masked_data[mask]["param"]}_mask'])

    DLA_mean_gz=masked_data['DLA']['data'][planes[plane]['index']].compressed()
    subDLA_mean_gz=masked_data['subDLA']['data'][planes[plane]['index']].compressed()
    LymanLimit_mean_gz=masked_data['LymanLimit']['data'][planes[plane]['index']].compressed()

    scatter=ax_sightlinescatter.scatter(nH_col,mean_gz,marker='x',s=2,alpha=.2,zorder=12,c=bin_radii,cmap='plasma_r')

    scatter_colourbar=fig_sightlinescatter.colorbar(scatter, ax=ax_sightlinescatter,location='left')
    scatter_colourbar.solids.set_alpha(1)
    scatter_colourbar.ax.tick_params(labelsize=16)
    scatter_colourbar.set_label('$R_{200_{crit}}$-Normalised Radial Distance From Centre Of FoF Group ($R_{200_{crit}}$)',fontsize=20)

    plt.xlabel(f'Pixel Projected ${plane}$ Planar Number Density ($H_1^1/cm^2$)',fontsize=18)
    plt.ylabel('Solar-Relative Pixel-Mass-Weighted Mean Metallicity ($Z_\odot$)',fontsize=18)
    ax_sightlinescatter.tick_params(labelsize=16)

    xlims=[np.float64(10**9),np.float64(10**23)]
    ylims=[np.float64(10**-6),np.float64(10**2)]
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(xlims)
    plt.ylim(ylims)

    ax_sightlinescatter.axvline(np.float64(10**20.3),c='r',ls='dashed')
    ax_sightlinescatter.fill_betweenx(np.array([10**-8,10**4]),np.float64(10**20.3),np.float64(10**23),color='r',alpha=.2)
    ax_sightlinescatter.text(np.float64(0.2*10**22),np.float64(10**-2.5),'DLA',c='r',rotation=45,fontsize=14)

    ax_sightlinescatter.axvline(10**19,c='b',ls='dashed')
    ax_sightlinescatter.fill_betweenx(np.array([10**-8,10**4]),np.float64(10**19),np.float64(10**20.3),color='b',alpha=.2)
    ax_sightlinescatter.text(np.float64(1.3*10**19),np.float64(10**-2.5),'Sub-DLA',c='b',rotation=45,fontsize=14)
    
    ax_sightlinescatter.axvline(10**17.2,c='g',ls='dashed')
    ax_sightlinescatter.fill_betweenx(np.array([10**-8,10**4]),np.float64(10**17.2),np.float64(10**19),color='g',alpha=.2)
    ax_sightlinescatter.text(np.float64(2.1*10**17),np.float64(10**-2.5),'Lyman Limit',c='g',rotation=45,fontsize=14)

    ax_sightlinescatter.axhline(10**-3,color='blueviolet',ls='dashed',lw=2)
    ax_sightlinescatter.text(np.float64(0.8*10**14),np.float64(1.2*10**-3),'Low Metallicity Threshold, $Z \leq 10^{-3}Z_{\odot}$',fontsize=18,color='blueviolet')

    ax_sightlinescatter.axhline(1,color='indigo',ls='dashed',lw=2)
    ax_sightlinescatter.text(np.float64(0.8*10**18),np.float64(1.2),'Solar Metallicity, $Z \leq Z_{\odot}$',fontsize=18,color='indigo')

    left,bottom,width,height=ax_sightlinescatter.get_position().bounds

    axis_bins=100

    ax_projdenshist=fig_sightlinescatter.add_axes([left,bottom+height,width,height/4],sharex=ax_sightlinescatter)
    
    plt.hist(nH_col,bins=np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),axis_bins),zorder=13,color='white',edgecolor='black',log=True)
    
    plt.yscale('log')
    plt.ylim([0,np.float64(10**6)])
    plt.title(f'{display_halo} {plane} Projection, z =${display_redshift}$, {bin_num} bins',fontsize=20,pad=20)

    ax_projdenshist.tick_params(labelbottom=False,labelleft=False,labelright=True,labelsize=16)
    ax_projdenshist.yaxis.tick_right()
    ax_projdenshist.set_ylabel('Number Density',fontsize=18)

    ax_projdenshist.axvline(np.float64(10**20.3),c='r',ls='dashed',zorder=11)
    ax_projdenshist.fill_betweenx([np.float64(0),np.float64(10**6)],np.float64(10**20.3),np.float64(10**23),color='r',alpha=.2)

    ax_projdenshist.axvline(10**19,c='b',ls='dashed',zorder=11)
    ax_projdenshist.fill_betweenx([np.float64(0),np.float64(10**6)],np.float64(10**19),np.float64(10**20.3),color='b',alpha=.2)

    ax_projdenshist.axvline(10**17.2,c='g',ls='dashed',zorder=11)
    ax_projdenshist.fill_betweenx([np.float64(0),np.float64(10**6)],np.float64(10**17.2),np.float64(10**19),color='g',alpha=.2,zorder=1)


    ax_gzhist=fig_sightlinescatter.add_axes([left+width,bottom,width/4,height],sharey=ax_sightlinescatter)
    plt.hist(mean_gz,bins=np.logspace(np.log10(ylims[0]),np.log10(ylims[1]),axis_bins),orientation='horizontal',color='white',edgecolor='black',log=True,zorder=13)
    
    plt.hist(subDLA_mean_gz,bins=np.logspace(np.log10(ylims[0]),np.log10(ylims[1]),axis_bins),orientation='horizontal',color='blue',edgecolor='cyan',alpha=.6,log=True,lw=2,zorder=13)
    plt.hist(LymanLimit_mean_gz,bins=np.logspace(np.log10(ylims[0]),np.log10(ylims[1]),axis_bins),orientation='horizontal',color='green',edgecolor='springgreen',alpha=.5,log=True,lw=2,zorder=13)
    plt.hist(DLA_mean_gz,bins=np.logspace(np.log10(ylims[0]),np.log10(ylims[1]),axis_bins),orientation='horizontal',color='red',edgecolor='orangered',alpha=.4,log=True,lw=2,zorder=13)

    ax_gzhist.axhline(10**-3,color='blueviolet',ls='dashed',lw=2)
    ax_gzhist.axhline(1,color='indigo',ls='dashed',lw=2)
    ax_gzhist.set_xlim([0,np.float64(10**6)])
    ax_gzhist.set_xscale('log')

    ax_gzhist.tick_params(labelleft=False,labelbottom=False,labeltop=True,labelsize=16)
    ax_gzhist.xaxis.tick_top()
    ax_gzhist.set_xlabel('Number Density',fontsize=18)
    gzhist_labels=ax_gzhist.get_xticklabels()
    gzhist_labels[0].set_visible(False)

    radii=['half','1','2','5']
    colours=['darkviolet','blue','dodgerblue','springgreen']

    for radius in radii:
        rad_nH_col=np.load(f'halos/{halo}/sightlines/{snap_num}/{radius}/nH_col.npy').flatten()
        rad_gz=np.load(f'halos/{halo}/sightlines/{snap_num}/{radius}/mean_gz.npy').flatten()
        ax_sightlinescatter.scatter(rad_nH_col,rad_gz,marker='x',s=2,alpha=.2,zorder=10,c=colours[radii.index(radius)])

        ax_projdenshist.hist(rad_nH_col,bins=np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),axis_bins),zorder=10,color=colours[radii.index(radius)],edgecolor=colours[radii.index(radius)],alpha=.4,log=True)
        ax_gzhist.hist(rad_gz,bins=np.logspace(np.log10(ylims[0]),np.log10(ylims[1]),axis_bins),orientation='horizontal',color=colours[radii.index(radius)],edgecolor=colours[radii.index(radius)],alpha=.4,log=True)

    if os.path.isdir(f'figures/{halo}/sightlines/{snap_num}') != True:
        os.makedirs(f'figures/{halo}/sightlines/{snap_num}',exist_ok=True)

    plt.savefig(f'figures/{halo}/sightlines/{snap_num}/sightline_scatters.png',format="png",dpi=250,bbox_inches='tight')



    plt.show()

def sightline_hists(halo,plane,**kwargs):
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

    snap_num=str(snap_num)
    while len(snap_num)<3: #Reformats snapshot number correctly into 3 digit string
        snap_num='0'+snap_num
        
    bin_num=512

    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    masked_data={'DLA':{'param':'mean_gz'},'subDLA':{'param':'mean_gz'},'LymanLimit':{'param':'mean_gz'}}

    for mask in masked_data:
        loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked/{mask}.npz')
        
        masked_data[mask]['obj'] = np.ma.masked_array(loaded[f'{masked_data[mask]["param"]}_data'], mask=loaded[f'{masked_data[mask]["param"]}_mask'])

    radii=['half','1','2','5']

    lims=[np.float64(10**-2),np.float64(10**2)]
    axis_bins=100

    for radius in radii:
        for mask in masked_data:
            loaded = np.load(f'halos/{halo}/sightlines/{snap_num}/{radius}/{mask}.npz')
            masked_data[mask][radius] = np.ma.masked_array(loaded[f'{masked_data[mask]["param"]}_data'], mask=loaded[f'{masked_data[mask]["param"]}_mask'])
        
    fig_sightlinehists,ax_sightlinehists=plt.subplots(3,figsize=(15,15),constrained_layout=True)

    row_index=0
    

    colours=['white','darkviolet','blue','dodgerblue','springgreen']
    row_colours=['red','blue','green']

    for mask in masked_data:
        colour_index=0
        ax_sightlinehists[row_index].set_title(mask)
        for radius in masked_data[mask]:
            if radius!='param':
                ax_sightlinehists[row_index].hist(masked_data[mask][radius][planes[plane]['index']].compressed(),bins=np.logspace(np.log10(lims[0]),np.log10(lims[1]),axis_bins),edgecolor=row_colours[row_index],color=colours[colour_index],alpha=.6,log=True,lw=2)
                colour_index+=1
        row_index+=1

    for ax in ax_sightlinehists:
        ax.set_xscale('log')
        ax.set_ylabel('Number Density',fontsize=18)
        ax.set_xlabel('Solar-Relative Pixel-Mass-Weighted Mean Metallicity ($Z_\odot$)',fontsize=18)

    if os.path.isdir(f'figures/{halo}/sightlines/{snap_num}') != True:
        os.makedirs(f'figures/{halo}/sightlines/{snap_num}',exist_ok=True)

    plt.savefig(f'figures/{halo}/sightlines/{snap_num}/sightline_hists.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()


