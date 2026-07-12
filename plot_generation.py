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

    type_plot_info={'gas':{'cmap':'plasma','title':f'Gas'},'dm':{'cmap':'viridis','title':f'Dark Matter'},'stars':{'cmap':'magma','title':f'Stars'}}

    fig_masshist,ax_masshist=plt.subplots(len(req_types),3,figsize=(15,13),constrained_layout=True)
    
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
            ax_masshist[row_index][planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(8.5, 16),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))
            ax_masshist[row_index][planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(16, 8.5),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))                

            ax_masshist[row_index][planes[plane]['index']].text(12.25,9.5,planes[plane]['x_label'].split()[0],fontsize=18,ha='center')
            ax_masshist[row_index][planes[plane]['index']].text(9.5,12,planes[plane]['y_label'].split()[0],fontsize=18,rotation='vertical',ha='center')

            ax_masshist[row_index][planes[plane]['index']].text(11.75,6.5,planes[plane]['x_label'].split()[1],fontsize=18,ha='center')
            ax_masshist[row_index][planes[plane]['index']].text(7,10.25,planes[plane]['y_label'].split()[1],fontsize=18,rotation='vertical',ha='center')

            if planes[plane]['index']==0:
                ax_masshist[row_index][planes[plane]['index']].set_ylabel(type_plot_info[matter_type]['title'],fontsize=24)
        
        colourbar=fig_masshist.colorbar(imshows[-1],ax=ax_masshist[row_index],shrink=1)
        colourbar.set_label('$log_{10}($Projected Density$)$ ($g/cm^2$)',fontsize=18)
        colourbar.ax.tick_params(labelsize=16)
        
        for ax in ax_masshist[row_index]:
            #ax.set_aspect('equal', adjustable='box')
            
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(18)
            #ax.yaxis.label.set_size(18)

            ax.set_xlim(-17,17)
            ax.set_ylim(-17,17)

            ax.tick_params(labelsize=16, axis='both',direction='in',top=True,right=True)
            ax.tick_params(axis='x', direction='in', pad=-20)
            ax.tick_params(axis='y', direction='in', pad=-30)


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
        if matter_type=='dm':
            def logline(rho0,r,n):
                return(rho0*r**n)
            print(densities[0])
            ax_radialdensity.plot(bin_centres[bin_centres<5],logline(3*10**-24,bin_centres[bin_centres<5],-1),ls='dashed',c='b',alpha=0.5)
            ax_radialdensity.plot(bin_centres[bin_centres>0.5],logline(6*10**-24,bin_centres[bin_centres>0.5],-3),ls='dashed',c='b',alpha=0.5)
            plt.text(2,logline(3*10**-24,4,-1)+0.1*10**-24,r'$\rho_{DM} \propto r^{-1}$',fontsize=12,c='b',rotation=-15)
            plt.text(20,logline(6*10**-24,50,-3)+0.2*10**-28,r'$\rho_{DM} \propto r^{-3}$',fontsize=12,c='b',rotation=-38)
    
    ax_radialdensity.axvline(halo_r200,c='k',ls='dashed')
    plt.text(halo_r200+1, 10**-26, '$R_{200_{crit}}=$'+str(np.round(halo_r200,1))+'$ kpc$', fontsize=14,c='k',rotation='vertical')

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Radius ($kpc$)',fontsize=14)
    plt.ylabel('Spherical Radial Density ($g/cm^3$)',fontsize=14)
    plt.xlim(left=0.05)
    plt.ylim(top=10**-22)
    #plt.title(display_halo+', $z=$'+str(display_redshift),fontsize=16)
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
        req_dens=['nH_col']

    dens_plot_info={'total_mass':{'cmap':'plasma','title':f'{display_halo} Total Gas, z =${display_redshift}$, {bin_num} bins','cbar_label':'$log_{10}($Projected Density$)$ ($g/cm^2$)'},'hydrogen_mass':{'cmap':'viridis','title':f'{display_halo} Hydrogen Mass Fraction, z =${display_redshift}$, {bin_num} bins','cbar_label':'$log_{10}($Projected Density$)$ ($g/cm^2$)'},'nH_col':{'cmap':'cividis','title':f'{display_halo} Neutral Hydrogen Fraction, z =${display_redshift}$, {bin_num} bins','cbar_label':'$log_{10}($Column Density$)$ ($atoms/cm^2$)'}}

    fig_gashist,ax_gashist=plt.subplots(1,3,figsize=(15,5),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    for density in req_dens:
        loaded_data={f'{plane}':np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/{density}.npy')[planes[plane]['index']]for plane in planes}
        loaded_data['rel_pos']=read_raw_file(halo,'gas','rel_pos',snap_num=snap_num)
        
        vmin=np.log10(np.min([np.min(loaded_data[plane][loaded_data[plane]!=0.0]) for plane in planes])) #Find minimum (excluding zeros) and maximum projected density values across projection axes to normalise colour bars to
        vmax=np.log10(np.max([loaded_data[plane]for plane in planes]))

        extents=calc.get_extent(loaded_data['rel_pos'])
        plot_extents={plane:[extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc), extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)] for plane in planes}

        row_index=req_dens.index(density)

        imshows=[ax_gashist[planes[plane]['index']].imshow(np.log10(loaded_data[plane]),vmin=vmin,vmax=vmax,extent=plot_extents[plane],cmap=dens_plot_info[density]['cmap'],aspect='equal') for plane in planes]
        for plane in planes:
            ax_gashist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(8.5, 16),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))
            ax_gashist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(16, 8.5),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))                

            ax_gashist[planes[plane]['index']].text(12.25,9.5,planes[plane]['x_label'].split()[0],fontsize=18,ha='center')
            ax_gashist[planes[plane]['index']].text(9.5,12,planes[plane]['y_label'].split()[0],fontsize=18,rotation='vertical',ha='center')

            ax_gashist[planes[plane]['index']].text(11.75,6.5,planes[plane]['x_label'].split()[1],fontsize=18,ha='center')
            ax_gashist[planes[plane]['index']].text(7,10.25,planes[plane]['y_label'].split()[1],fontsize=18,rotation='vertical',ha='center')
        
        colourbar=fig_gashist.colorbar(imshows[-1],ax=ax_gashist,shrink=.9)

        colourbar.ax.axhline(20.3,c='r',lw=2)
        colourbar.ax.axhline(19,c='b',lw=2)
        colourbar.ax.axhline(17.2,c='g',lw=2)

        colourbar.ax.fill_betweenx([20.3,24],0,1,facecolor='none',edgecolor='red',hatch='///')
        colourbar.ax.fill_betweenx([19,20.3],0,1,facecolor='none',edgecolor='blue',hatch='///')
        colourbar.ax.fill_betweenx([17.2,19],0,1,facecolor='none',edgecolor='green',hatch='///')


        '''colourbar.ax.text(-2.2,20.3,'DLA',fontsize=16,c='r',rotation=45)
        colourbar.ax.text(-3.7,18,'Sub-DLA',fontsize=16,c='b',rotation=45)
        colourbar.ax.text(-2,17.5,'LLS',fontsize=16,c='g',rotation=45)'''
        colourbar.set_label(dens_plot_info[density]['cbar_label'],fontsize=16,ha='center')
        colourbar.ax.tick_params(labelsize=16)
        
        for ax in ax_gashist:
            #ax.set_aspect('equal', adjustable='box')
            
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(18)
            #ax.yaxis.label.set_size(18)

            ax.set_xlim(-17,17)
            ax.set_ylim(-17,17)

            ax.tick_params(labelsize=16, axis='both',direction='in',top=True,right=True)
            ax.tick_params(axis='x', direction='in', pad=-20)
            ax.tick_params(axis='y', direction='in', pad=-30)

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

    fig_gzhist,ax_gzhist=plt.subplots(1,3,figsize=(15,5),constrained_layout=True)
    
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
            ax_gzhist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(8.5, 16),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))
            ax_gzhist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(16, 8.5),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))                

            ax_gzhist[planes[plane]['index']].text(12.25,9.5,planes[plane]['x_label'].split()[0],fontsize=18,ha='center')
            ax_gzhist[planes[plane]['index']].text(9.5,12,planes[plane]['y_label'].split()[0],fontsize=18,rotation='vertical',ha='center')

            ax_gzhist[planes[plane]['index']].text(11.75,6.5,planes[plane]['x_label'].split()[1],fontsize=18,ha='center')
            ax_gzhist[planes[plane]['index']].text(7,10.25,planes[plane]['y_label'].split()[1],fontsize=18,rotation='vertical',ha='center')
        
    colourbar=fig_gzhist.colorbar(imshows[-1],ax=ax_gzhist,shrink=.9)
    colourbar.set_label('$log_{10}($Mass-Weighted Mean Solar-Relative Metallicity$)$ ($Z_\odot$)',fontsize=12)
    colourbar.ax.tick_params(labelsize=16)
        
    for ax in ax_gzhist:
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(18)
            #ax.yaxis.label.set_size(18)

            ax.set_xlim(-17,17)
            ax.set_ylim(-17,17)

            ax.tick_params(labelsize=16, axis='both',direction='in',top=True,right=True)
            ax.tick_params(axis='x', direction='in', pad=-20)
            ax.tick_params(axis='y', direction='in', pad=-30)

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
    ax_scatter.text(np.float64(6*10**17),np.float64(10**-2.5),'LLS',c='g',rotation=45,fontsize=14)

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
        #plt.title(f'{display_halo} {plane} Projection, z =${display_redshift}$, {bin_num} bins',fontsize=20,pad=20)

        ax_projdenshist.tick_params(labelbottom=False,labelleft=False,labelright=True,labelsize=16)
        ax_projdenshist.yaxis.tick_right()
        ax_projdenshist.set_ylabel('$log_{10}(N_{px})$',fontsize=18)

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
        ax_gzhist.set_xlim([0,np.float64(10**5)])
        ax_gzhist.set_xscale('log')

        ax_gzhist.tick_params(labelleft=False,labelbottom=False,labeltop=True,labelsize=16)
        ax_gzhist.xaxis.tick_top()
        ax_gzhist.set_xlabel('$log_{10}(N_{px})$',fontsize=18)
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
    #plt.title(f'{display_halo}, {plane}-plane, {bin_num} bins',fontsize=20)
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

    fig_gashist,ax_gashist=plt.subplots(1,3,figsize=(15,5),constrained_layout=True)
    
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
            ax_gashist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(8.5, 16),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))
            ax_gashist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(16, 8.5),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))                

            ax_gashist[planes[plane]['index']].text(12.25,9.5,planes[plane]['x_label'].split()[0],fontsize=18,ha='center')
            ax_gashist[planes[plane]['index']].text(9.5,12,planes[plane]['y_label'].split()[0],fontsize=18,rotation='vertical',ha='center')

            ax_gashist[planes[plane]['index']].text(11.75,6.5,planes[plane]['x_label'].split()[1],fontsize=18,ha='center')
            ax_gashist[planes[plane]['index']].text(7,10.25,planes[plane]['y_label'].split()[1],fontsize=18,rotation='vertical',ha='center')
    
    colourbar=fig_gashist.colorbar(imshows[-1],ax=ax_gashist,shrink=.9)
    colourbar.set_label('$log_{10}($Projected Density$)$ ($g/cm^2$)',fontsize=18)
    colourbar.ax.tick_params(labelsize=16)
    
    for ax in ax_gashist:
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(18)
            #ax.yaxis.label.set_size(18)

            ax.set_xlim(-17,17)
            ax.set_ylim(-17,17)

            ax.tick_params(labelsize=16, axis='both',direction='in',top=True,right=True)
            ax.tick_params(axis='x', direction='in', pad=-20)
            ax.tick_params(axis='y', direction='in', pad=-30)

    if os.path.isdir(f'figures/{halo}/{bin_num}/contour_gas_densities')!=True:
        os.makedirs(f'figures/{halo}/{bin_num}/contour_gas_densities',exist_ok=True)

    plt.savefig(f'figures/{halo}/{bin_num}/contour_gas_densities/{snap_num}.pdf',format="pdf",dpi=250,bbox_inches='tight')
    
    plt.show()


def stellar_masses_redshift(halos):
    if os.path.isdir(f'figures/all_halos/stellar_masses')!=True:
        os.makedirs(f'figures/all_halos/stellar_masses')

    fig_masses,ax_masses=plt.subplots(figsize=(16,6))

    colours=['lightseagreen','deeppink','darkorange']
    markers=['D','D','D']

    nH_snaps={'T1_Aug':{'DLA':[[0,5],[6,7]],'subDLA':[[0,9],[11,13]],'LymanLimit':[[0,13]]},'halo8':{'DLA':[[0,1],[4,19],[26,29],[35,38]],'subDLA':[[0,38]],'LymanLimit':[[0,38]]},'T4_Aug':{'DLA':[[1,2],[3,4]],'subDLA':[[0,12]],'LymanLimit':[[0,13]]}}

    formation_events={'T1_Aug':{'Smaller Merger Begins':(2,'dashed'),'Small Merger Concludes, Main Merger Begins':(3,'dashed'),'Gas Merger Concludes':(5,'dashed'),'Overall Merger Concludes':(7,'dashed')},'T4_Aug':{'Merger Begins':(3,'dashed'),'Merger Concludes':(6,'dashed')},'halo8':{'Satellite In Halo':(15,'dotted'),'Central Disturbance':(19,'dotted'),'Disturbance Resolved':(21,'dotted'),'1st Merger Begins':(22,'dashed'),'1st Merger Concluding':(28,'dashed'),'2nd Merger In Progress':(31,'dashed'),'2nd Gas Merger Concludes':(32,'dashed'),'Tertiary Object In Halo':(33,'dotted'),'2nd Merger Concludes, 3rd Merger Begins':(34,'dashed'),'3rd Merger Concluding':(37,'dashed')}}


    for halo in halos:
        display_halo=halo.replace('_',' ')

        redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
        masses=np.load(f'halos/{halo}/stellar_masses/stellar_masses.npy')
          
        ax_masses.plot(redshifts,masses,c=colours[halos.index(halo)],zorder=1,lw=2)

        DLA_snaps=nH_snaps[halo]['DLA']
        subDLA_snaps=nH_snaps[halo]['subDLA']
        LLS_snaps=nH_snaps[halo]['LymanLimit']

        DLA_points=[ax_masses.scatter(redshifts[snaps[0]:snaps[1]],masses[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_masses.scatter(redshifts[snaps[0]:snaps[1]],masses[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_masses.scatter(redshifts[snaps[0]:snaps[1]],masses[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]

        events=[ax_masses.axvline(redshifts[formation_events[halo][event][0]],c=colours[halos.index(halo)],linestyle=formation_events[halo][event][1],alpha=0.5,zorder=0) for event in formation_events[halo]]

        event_labels=[ax_masses.text(redshifts[formation_events[halo][event][0]]+0.03,1.4*10**7,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor') if event=='Merger Begins' else ax_masses.text(redshifts[formation_events[halo][event][0]],1.3*10**7,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor')for event in formation_events[halo]]

        if halo=='halo8':
            shift=0.15
        else:
            shift=0.2
        plt.text(redshifts[0]+shift,masses[0],display_halo,fontsize=14,c=colours[halos.index(halo)])
    

    ax_masses.invert_xaxis()

    ax_masses.xaxis.label.set_size(18)
    ax_masses.yaxis.label.set_size(18)
    ax_masses.tick_params(labelsize=18)
    plt.xlabel('Redshift (z)',fontsize=18)
    plt.ylabel(r'Total Stellar Mass ($M_{\odot}$)',fontsize=18)
    plt.xlim([4.3,0.9])
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

    fig_zs,ax_zs=plt.subplots(figsize=(16,6))

    colours=['lightseagreen','deeppink','darkorange']
    markers=['D','D','D']

    nH_snaps={'T1_Aug':{'DLA':[[0,5],[6,7]],'subDLA':[[0,9],[11,13]],'LymanLimit':[[0,13]]},'halo8':{'DLA':[[0,1],[4,19],[26,29],[35,38]],'subDLA':[[0,38]],'LymanLimit':[[0,38]]},'T4_Aug':{'DLA':[[1,2],[3,4]],'subDLA':[[0,12]],'LymanLimit':[[0,13]]}}

    formation_events={'T1_Aug':{'Smaller Merger Begins':(2,'dashed'),'Small Merger Concludes, Main Merger Begins':(3,'dashed'),'Gas Merger Concludes':(5,'dashed'),'Overall Merger Concludes':(7,'dashed')},'T4_Aug':{'Merger Begins':(3,'dashed'),'Merger Concludes':(6,'dashed')},'halo8':{'Satellite In Halo':(15,'dotted'),'Central Disturbance':(19,'dotted'),'Disturbance Resolved':(21,'dotted'),'1st Merger Begins':(22,'dashed'),'1st Merger Concluding':(28,'dashed'),'2nd Merger In Progress':(31,'dashed'),'2nd Gas Merger Concludes':(32,'dashed'),'Tertiary Object In Halo':(33,'dotted'),'2nd Merger Concludes, 3rd Merger Begins':(34,'dashed'),'3rd Merger Concluding':(37,'dashed')}}

    for halo in halos:
        display_halo=halo.replace('_',' ')

        redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
        gas_zs=np.load(f'halos/{halo}/metallicity_behaviour/gas.npy')
        star_zs=np.load(f'halos/{halo}/metallicity_behaviour/stars.npy')

        if halo == 'T4_Aug':
            snap_0_gas_z=(0.0033585230966704242/z_sol+gas_zs[1])/2
            gas_zs[0]=snap_0_gas_z
          
        ax_zs.plot(redshifts,gas_zs,c=colours[halos.index(halo)],zorder=1,lw=2)
        ax_zs.plot(redshifts,star_zs,c=colours[halos.index(halo)],linestyle='dashdot',zorder=0,lw=2)

        print(gas_zs)

        DLA_snaps=nH_snaps[halo]['DLA']
        subDLA_snaps=nH_snaps[halo]['subDLA']
        LLS_snaps=nH_snaps[halo]['LymanLimit']

        DLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],gas_zs[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],gas_zs[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],gas_zs[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]

        DLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],star_zs[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],star_zs[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],star_zs[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]

        events=[ax_zs.axvline(redshifts[formation_events[halo][event][0]],c=colours[halos.index(halo)],linestyle=formation_events[halo][event][1],alpha=0.5,zorder=0) for event in formation_events[halo]]

        event_labels=[ax_zs.text(redshifts[formation_events[halo][event][0]]+0.03,0.245,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor') if event=='Merger Begins' else ax_zs.text(redshifts[formation_events[halo][event][0]],0.237,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor')for event in formation_events[halo]]

        if halo=='halo8':
            shift=0.31
        else:
            shift=0.35
        plt.text(redshifts[0]-0.03+shift,gas_zs[0],f'{display_halo} Gas',fontsize=14,c=colours[halos.index(halo)])
        plt.text(redshifts[0]+shift,star_zs[0],f'{display_halo} Stars',fontsize=14,c=colours[halos.index(halo)])

    ax_zs.invert_xaxis()

    ax_zs.xaxis.label.set_size(18)
    ax_zs.yaxis.label.set_size(18)
    ax_zs.tick_params(labelsize=18)
    plt.xlabel('Redshift (z)',fontsize=18)
    plt.ylabel(r'Metallicity ($Z_\odot$)',fontsize=18)
    plt.xlim([4.4,0.9])
    plt.yscale('log')

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

    fig_zs,ax_zs=plt.subplots(figsize=(16,6))

    colours=['lightseagreen','deeppink','darkorange']
    markers=['D','D','D']

    nH_snaps={'T1_Aug':{'DLA':[[0,5],[6,7]],'subDLA':[[0,9],[11,13]],'LymanLimit':[[0,13]]},'halo8':{'DLA':[[0,1],[4,19],[26,29],[35,38]],'subDLA':[[0,38]],'LymanLimit':[[0,38]]},'T4_Aug':{'DLA':[[1,2],[3,4]],'subDLA':[[0,12]],'LymanLimit':[[0,13]]}}

    formation_events={'T1_Aug':{'Smaller Merger Begins':(2,'dashed'),'Small Merger Concludes, Main Merger Begins':(3,'dashed'),'Gas Merger Concludes':(5,'dashed'),'Overall Merger Concludes':(7,'dashed')},'T4_Aug':{'Merger Begins':(3,'dashed'),'Merger Concludes':(6,'dashed')},'halo8':{'Satellite In Halo':(15,'dotted'),'Central Disturbance':(19,'dotted'),'Disturbance Resolved':(21,'dotted'),'1st Merger Begins':(22,'dashed'),'1st Merger Concluding':(28,'dashed'),'2nd Merger In Progress':(31,'dashed'),'2nd Gas Merger Concludes':(32,'dashed'),'Tertiary Object In Halo':(33,'dotted'),'2nd Merger Concludes, 3rd Merger Begins':(34,'dashed'),'3rd Merger Concluding':(37,'dashed')}}

    for halo in halos:
        display_halo=halo.replace('_',' ')

        redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
        m200s=np.load(f'halos/{halo}/m200s/m200.npy')

        DLA_snaps=nH_snaps[halo]['DLA']
        subDLA_snaps=nH_snaps[halo]['subDLA']
        LLS_snaps=nH_snaps[halo]['LymanLimit']

        DLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],m200s[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],m200s[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],m200s[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]

        events=[ax_zs.axvline(redshifts[formation_events[halo][event][0]],c=colours[halos.index(halo)],linestyle=formation_events[halo][event][1],alpha=0.5,zorder=0) for event in formation_events[halo]]

        event_labels=[ax_zs.text(redshifts[formation_events[halo][event][0]]+0.03,4.2*10**9,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor') if event=='Merger Begins' else ax_zs.text(redshifts[formation_events[halo][event][0]],4.1*10**9,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor')for event in formation_events[halo]]

        ax_zs.plot(redshifts,m200s,c=colours[halos.index(halo)],zorder=1,lw=2)
        if halo=='halo8':
            shift=0.15
        else:
            shift=0.2
        plt.text(redshifts[0]+shift,m200s[0],f'{display_halo}',fontsize=14,c=colours[halos.index(halo)])


    ax_zs.invert_xaxis()

    ax_zs.xaxis.label.set_size(18)
    ax_zs.yaxis.label.set_size(18)
    ax_zs.tick_params(labelsize=18)
    plt.xlabel('Redshift (z)',fontsize=18)
    plt.ylabel(r'$M_{200_{crit}} (M_\odot)$',fontsize=18)
    plt.xlim([4.3,0.9])
    plt.yscale('log')


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

    fig_weighted,ax_weighted=plt.subplots(figsize=(8,6))

    xlims=[np.float64(0.3*10**-3),np.float64(1)]

    plt.hist(unweighted_gz,edgecolor='blue',log=True,bins=np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),100),histtype='step')
    plt.hist(weighted_gz,color='dodgerblue',edgecolor='darkblue',log=True,bins=np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),100),alpha=.4)
    
    plt.text(0.09,150,'Unweighted Mean',c='b',fontsize=16)
    

    plt.xlim(xlims)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Pixel Metallicity Value ($Z_\odot$)',fontsize=16)
    plt.ylabel('$log_{10}(N_{px})$',fontsize=16)

    ax_weighted.xaxis.label.set_size(18)
    ax_weighted.yaxis.label.set_size(18)
    ax_weighted.tick_params(labelsize=18)

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

    nh_masks=['DLA']

    fig_sightlineconts,ax_sightlineconts=plt.subplots(figsize=(7,7),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'}}

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
            
            raw_contours=[ax_sightlineconts.contour(rad_masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=raw_extents[plane],colors=colours[radii.index(radius)],zorder=3) for plane in planes if np.any(~rad_masks[col_dens]['data'][planes[plane]['index']].mask)]

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

            contours=[ax_sightlineconts.contour(rad_masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=extents[plane],colors=colours[radii.index(radius)],zorder=3,linewidths=5,alpha=0.7) for plane in planes if np.any(~rad_masks[col_dens]['data'][planes[plane]['index']].mask)]
            
        masks={'DLA':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.3],[1,1,1,0]]))},'subDLA':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.3],[1,1,1,0]]))},'LymanLimit':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.3],[1,1,1,0]]))}}

        loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/512px/gas_only/masked/{col_dens}.npz') 
        masks[col_dens]['data'] = np.ma.masked_array(loaded[f'nH_col_data'], mask=loaded[f'nH_col_mask'])

        lo_z_loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/512px/gas_only/masked/{col_dens}.npz')
        masks[col_dens]['lo_z_data'] = np.ma.masked_array(lo_z_loaded[f'nH_col_data'], mask=lo_z_loaded[f'nH_col_mask'])

        raw_obj_contours=[ax_sightlineconts.contour(masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=obj_plot_extents[plane],colors=masks[col_dens]['colour']) for plane in planes if np.any(~masks[col_dens]['data'][planes[plane]['index']].mask)]

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

        adj_obj_contours=[ax_sightlineconts.contour(masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=shifted_obj_extents[plane],colors=masks[col_dens]['colour'],zorder=2,linewidths=5) for plane in planes if np.any(~masks[col_dens]['data'][planes[plane]['index']].mask)]
        adj_obj_fills=[ax_sightlineconts.imshow(masks[col_dens]['data'][planes[plane]['index']].mask.astype(float),extent=shifted_obj_extents[plane],cmap=masks[col_dens]['cmap'],zorder=1) for plane in planes if np.any(~masks[col_dens]['data'][planes[plane]['index']].mask)]  
        adj_obj_imshows=[ax_sightlineconts.imshow(np.log10(loaded_data[plane]),extent=shifted_obj_extents[plane],vmin=vmin,vmax=vmax,cmap=dens_plot_info['cmap'],aspect='equal',zorder=0) for plane in planes]
       
        for plane in planes:
            ax_sightlineconts.set_xlabel(planes[plane]['x_label'])
            ax_sightlineconts.set_ylabel(planes[plane]['y_label'])
            if planes[plane]['index']==1:
                ax_sightlineconts.set_title(col_dens,fontsize=20,pad=20)
    
    radii=[0.5,1,2,5,10]

    row_dims=[1,1.5,2.5]
    row_index=0

    ax_sightlineconts.set_xlim(-row_dims[row_index],row_dims[row_index])
    ax_sightlineconts.set_ylim(-row_dims[row_index],row_dims[row_index])
    ax_sightlineconts.set_aspect('equal', adjustable='box')
    ax_sightlineconts.set_aspect('equal', adjustable='box')
    ax_sightlineconts.set_box_aspect(1.0)
    ax_sightlineconts.xaxis.label.set_size(18)
    ax_sightlineconts.yaxis.label.set_size(18)

    ax_sightlineconts.yaxis.set_label_coords(-0.1, 0.5)

    ax_sightlineconts.tick_params(labelsize=16)

    ax_radscale=fig_sightlineconts.add_axes([0.965,0.1,0.12,0.833])


    ax_radscale.yaxis.tick_right()
    ax_radscale.set_ylim(0,5.5)
    ax_radscale.yaxis.set_label_position('right')
    ax_radscale.set_ylabel(r'Sample Annulus Radius ($R_{200_{crit}}$)',fontsize=14)
    ax_radscale.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax_radscale.xaxis.set_tick_params(labelsize=12)

    fractions=[0.5,1,2,5]
    labels=[r'$\frac{1}{2}R_{200_{crit}}$',r'$R_{200_{crit}}$',r'$2R_{200_{crit}}$',r'$5R_{200_{crit}}$']

    for frac in fractions:
        ax_radscale.axhline(frac,c=colours[fractions.index(frac)])
        ax_radscale.text(0.1,frac+0.1,labels[fractions.index(frac)],c=colours[fractions.index(frac)],fontsize=16)
    

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


def report_subfind_redshifts(halos):
    if os.path.isdir(f'figures/all_halos/subfinds')!=True:
        os.makedirs(f'figures/all_halos/subfinds')

    fig_subfinds, ax_subfinds=plt.subplots(3,figsize=(10,6),sharex=True)

    colours=['lightseagreen','deeppink','darkorange']
    markers=['D','D','D']

    nH_snaps={'T1_Aug':{'DLA':[[0,5],[6,7]],'subDLA':[[0,9],[11,13]],'LymanLimit':[[0,13]]},'halo8':{'DLA':[[0,1],[4,19],[26,29],[35,38]],'subDLA':[[0,38]],'LymanLimit':[[0,38]]},'T4_Aug':{'DLA':[[1,2],[3,4]],'subDLA':[[0,12]],'LymanLimit':[[0,13]]}}

    display_names=['T1','T4','h8']

    formation_events={'T1_Aug':{'$T1_A$':(2,'dashed'),'$T1_B$':(3,'dashed'),'$T1_C$':(4,'dotted'),'$T1_D$':(5,'dotted'),'$T1_E$':(7,'dashed')},'T4_Aug':{'$T4_A$':(3,'dashed'),'$T4_B$':(6,'dashed')},'halo8':{'$h8_A$':(18,'dotted'),'$h8_B$':(19,'dotted'),'$h8_C$':(22,'dashed'),'$h8_D$':(24,'dashed'),'$h8_E$':(29,'dashed'),'$h8_F$':(33,'dashed'),'$h8_G$':(35,'dashed')}}

    DLA_windows={'T1_Aug':{'start':0,'end':5},'T4_Aug':{'start':0,'end':0},'halo8':{'start':8,'end':19}}


    for halo in halos:

        redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
        m200s=np.load(f'halos/{halo}/m200s/m200.npy')

        DLA_snaps=nH_snaps[halo]['DLA']
        subDLA_snaps=nH_snaps[halo]['subDLA']
        LLS_snaps=nH_snaps[halo]['LymanLimit']

        '''DLA_points=[ax_subfinds[0].scatter(redshifts[snaps[0]:snaps[1]],m200s[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_subfinds[0].scatter(redshifts[snaps[0]:snaps[1]],m200s[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_subfinds[0].scatter(redshifts[snaps[0]:snaps[1]],m200s[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]'''

        events=[ax_subfinds[0].axvline(redshifts[formation_events[halo][event][0]],c=colours[halos.index(halo)],linestyle=formation_events[halo][event][1],alpha=0.5,zorder=0) for event in formation_events[halo]]

        event_labels=[ax_subfinds[0].text(redshifts[formation_events[halo][event][0]]-0.05,7.5*10**9,event,c=colours[halos.index(halo)],rotation='vertical') if event=='$T4_A$' else ax_subfinds[0].text(redshifts[formation_events[halo][event][0]]+0.02,7.5*10**9,event,c=colours[halos.index(halo)],rotation='vertical')for event in formation_events[halo]]

        ax_subfinds[0].plot(redshifts,m200s,c=colours[halos.index(halo)],zorder=1,lw=2)
        ax_subfinds[0].plot(redshifts[DLA_windows[halo]['start']:DLA_windows[halo]['end']],m200s[DLA_windows[halo]['start']:DLA_windows[halo]['end']],c='r',lw=8,alpha=0.3,zorder=0)
        
        if halo=='T4_Aug':
            shift=1
        else:
            shift=2 
        ax_subfinds[0].text(redshifts[-1]-0.02,m200s[-1]-shift*10**8,display_names[halos.index(halo)],color=colours[halos.index(halo)],fontsize=14)

        masses=np.load(f'halos/{halo}/stellar_masses/stellar_masses.npy')

        ax_subfinds[1].plot(redshifts,masses,c=colours[halos.index(halo)],zorder=1,lw=2)
        ax_subfinds[1].plot(redshifts[DLA_windows[halo]['start']:DLA_windows[halo]['end']],masses[DLA_windows[halo]['start']:DLA_windows[halo]['end']],c='r',lw=8,alpha=0.3,zorder=0)


        '''DLA_points=[ax_subfinds[1].scatter(redshifts[snaps[0]:snaps[1]],masses[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_subfinds[1].scatter(redshifts[snaps[0]:snaps[1]],masses[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_subfinds[1].scatter(redshifts[snaps[0]:snaps[1]],masses[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]'''

        events=[ax_subfinds[1].axvline(redshifts[formation_events[halo][event][0]],c=colours[halos.index(halo)],linestyle=formation_events[halo][event][1],alpha=0.5,zorder=0) for event in formation_events[halo]]

        if halo=='T1_Aug':
            shift=5
        else:
            shift=2 
        ax_subfinds[1].text(redshifts[-1]-0.02,masses[-1]-shift*10**5,display_names[halos.index(halo)],color=colours[halos.index(halo)],fontsize=14)

        
        gas_zs=np.load(f'halos/{halo}/metallicity_behaviour/gas.npy')
        star_zs=np.load(f'halos/{halo}/metallicity_behaviour/stars.npy')

        if halo == 'T4_Aug':
            snap_0_gas_z=(0.0033585230966704242/z_sol+gas_zs[1])/2
            gas_zs[0]=snap_0_gas_z
          
        ax_subfinds[2].plot(redshifts,gas_zs,c=colours[halos.index(halo)],zorder=1,lw=2)
        ax_subfinds[2].plot(redshifts[DLA_windows[halo]['start']:DLA_windows[halo]['end']],gas_zs[DLA_windows[halo]['start']:DLA_windows[halo]['end']],c='r',lw=8,alpha=0.3,zorder=0)
        ax_subfinds[2].plot(redshifts,star_zs,c=colours[halos.index(halo)],linestyle='dashdot',zorder=1,lw=2)
        ax_subfinds[2].plot(redshifts[DLA_windows[halo]['start']:DLA_windows[halo]['end']],star_zs[DLA_windows[halo]['start']:DLA_windows[halo]['end']],c='r',lw=8,alpha=0.3,zorder=0)

        '''DLA_points=[ax_subfinds[2].scatter(redshifts[snaps[0]:snaps[1]],gas_zs[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_subfinds[2].scatter(redshifts[snaps[0]:snaps[1]],gas_zs[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_subfinds[2].scatter(redshifts[snaps[0]:snaps[1]],gas_zs[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]

        DLA_points=[ax_subfinds[2].scatter(redshifts[snaps[0]:snaps[1]],star_zs[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_subfinds[2].scatter(redshifts[snaps[0]:snaps[1]],star_zs[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_subfinds[2].scatter(redshifts[snaps[0]:snaps[1]],star_zs[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]'''

        events=[ax_subfinds[2].axvline(redshifts[formation_events[halo][event][0]],c=colours[halos.index(halo)],linestyle=formation_events[halo][event][1],alpha=0.5,zorder=0) for event in formation_events[halo]]

        if halo=='T4_Aug':
            shift=0
            suffix=''
        elif halo=='halo8':
            shift=0.2
            suffix=' Gas'
        else:
            shift=2 
            suffix=''
        ax_subfinds[2].text(redshifts[-1]-0.02,gas_zs[-1]-shift*10**-2,display_names[halos.index(halo)]+suffix,color=colours[halos.index(halo)],fontsize=14)

        if halo=='halo8':
            shift=0.2
            ax_subfinds[2].text(redshifts[-1]-0.02,star_zs[-1]-shift*10**-2,display_names[halos.index(halo)]+' Stars',color=colours[halos.index(halo)],fontsize=14)


    ax_subfinds[0].set_ylim([ax_subfinds[0].get_ylim()[0],ax_subfinds[0].get_ylim()[1]+2*10**9])
    ax_subfinds[1].set_ylim([ax_subfinds[1].get_ylim()[0],ax_subfinds[1].get_ylim()[1]+0.5*10**7])


    for ax in ax_subfinds:
        ax.invert_xaxis()
        ax.set_xlim([4.1,0.6])
        ax.set_yscale('log')
        ax.tick_params(labelsize=12)



    
    ax_subfinds[0].tick_params(axis='x', top=False, bottom=False, labeltop=False, labelbottom=False, direction='in') 
    ax_subfinds[1].tick_params(axis='x', top=True, bottom=False, labeltop=False, labelbottom=False, direction='in') 
    ax_subfinds[2].tick_params(axis='x', bottom=True, labelbottom=True,direction='out')

    ax_subfinds[0].set_ylabel(r'$M_{200_{crit}} (M_\odot)$',fontsize=12)
    ax_subfinds[1].set_ylabel(r'Stellar Mass ($M_{\odot}$)',fontsize=12)
    ax_subfinds[2].set_ylabel(r'Metallicity ($Z_\odot$)',fontsize=12)
    ax_subfinds[2].set_xlabel('Redshift (z)',fontsize=12)

    ax_top = ax_subfinds[2].twiny()
    ax_top.set_xlim(ax_subfinds[2].get_xlim())

    ax_top.tick_params(axis='x', top=True, labeltop=False, direction='in')

    ax_top.spines['bottom'].set_visible(False)
    ax_top.spines['left'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    plt.subplots_adjust(hspace=0)
    
    plt.savefig('figures/all_halos/subfinds/subfinds.pdf',format="pdf",dpi=250,bbox_inches='tight')
    plt.show()

def report_T1_radial():
    
    snapshots=[152,155,156,157,159,164]
    colours=['red','orangered','gold','limegreen','dodgerblue','indigo']

    fig_radialdensity,ax_radialdensity=plt.subplots()

    halo='T1_Aug'
    matter_type='gas'

    for snap_num in snapshots:
        loaded_data=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/radial_mass_density/{matter_type}.npy')
        
        bin_centres=loaded_data[0]
        densities=loaded_data[1]

        ax_radialdensity.plot(bin_centres,densities,c=colours[snapshots.index(snap_num)],alpha=0.7)
        #plt.text(bin_centres[0], densities[0], snap_num, fontsize=14,c=colours[snapshots.index(snap_num)])
       
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Radius ($kpc$)',fontsize=14)
    plt.ylabel('Spherical Radial Density ($g/cm^3$)',fontsize=14)
    plt.xlim(left=0.1,right=50)
    plt.ylim(top=10**-22)
    #plt.title(display_halo+', $z=$'+str(display_redshift),fontsize=16)
    ax_radialdensity.xaxis.set_tick_params(labelsize=12)
    ax_radialdensity.yaxis.set_tick_params(labelsize=12)

    #ax_radialdensity.text(32,0.15*10**-34,'T1',c='k',fontsize=20)

    pos=ax_radialdensity.get_position()

    ax_radialdensity.set_position([pos.x0, pos.y0, pos.width, pos.height*0.9])

    ax_redshifts=fig_radialdensity.add_axes([pos.x0, pos.y0+pos.height*0.9, pos.width, pos.height*0.1])

    redshifts=[4.008,2.494,2.208,1.960,1.592,0.997]
    labels=[r'$z\approx4$',r'$T1_B$',r'$T1_C$',r'$T1_D$',r'$T1_E$',r'$z\approx1$']
    for redshift in redshifts:
        ax_redshifts.axvline(redshift,c=colours[redshifts.index(redshift)])
        ax_redshifts.text(redshift+0.15,-0.6,labels[redshifts.index(redshift)],c=colours[redshifts.index(redshift)],fontsize=14)
    
    ax_redshifts.invert_xaxis()
    ax_redshifts.xaxis.tick_top()
    ax_redshifts.xaxis.set_label_position('top')
    ax_redshifts.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    ax_redshifts.xaxis.set_tick_params(labelsize=12)
    ax_redshifts.set_xlabel('Redshift (z)',fontsize=14)
    ax_redshifts.set_xlim(4.2,0.8)


    plt.savefig('figures/T1_Aug/report_radials.pdf',format="pdf",dpi=250,bbox_inches='tight')


    plt.show()

def report_T4_densities():

    snapshots=[152,153]
    halo='T4_Aug'
    bin_num='512'
    matter_type='gas'
    titles=[f'$z=4.008$',f'$z=3.344$']

    type_plot_info={'gas':{'cmap':'plasma','title':f'Gas'},'dm':{'cmap':'viridis','title':f'Dark Matter'},'stars':{'cmap':'magma','title':f'Stars'}}

    fig_masshist,ax_masshist=plt.subplots(len(snapshots),3,figsize=(15,9),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    for snap_num in snapshots:
        loaded_data={f'{plane}':np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/total_mass/{matter_type}.npy')[planes[plane]['index']]for plane in planes}
        loaded_data['rel_pos']=read_raw_file(halo,matter_type,'rel_pos',snap_num=snap_num)
        
        vmin=np.log10(np.min([np.min(loaded_data[plane][loaded_data[plane]!=0.0]) for plane in planes])) #Find minimum (excluding zeros) and maximum projected density values across projection axes to normalise colour bars to
        vmax=np.log10(np.max([loaded_data[plane]for plane in planes]))

        extents=calc.get_extent(loaded_data['rel_pos'])
        plot_extents={plane:[extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc), extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)] for plane in planes}

        row_index=snapshots.index(snap_num)

        imshows=[ax_masshist[row_index][planes[plane]['index']].imshow(np.log10(loaded_data[plane]),extent=plot_extents[plane],vmin=vmin,vmax=vmax,cmap=type_plot_info[matter_type]['cmap'],aspect='equal') for plane in planes]
        for plane in planes:
            ax_masshist[row_index][planes[plane]['index']].annotate('', xytext=(5.5, 5.5), xy=(5.5, 11),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))
            ax_masshist[row_index][planes[plane]['index']].annotate('', xytext=(5.5, 5.5), xy=(11, 5.5),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))                

            ax_masshist[row_index][planes[plane]['index']].text(8.25,4.25,planes[plane]['x_label'].split()[0],fontsize=18,ha='center')
            ax_masshist[row_index][planes[plane]['index']].text(4.5,8,planes[plane]['y_label'].split()[0],fontsize=18,rotation='vertical',ha='center')

            ax_masshist[row_index][planes[plane]['index']].text(8.25,6.25,planes[plane]['x_label'].split()[1],fontsize=14,ha='center')
            ax_masshist[row_index][planes[plane]['index']].text(6.5,7.25,planes[plane]['y_label'].split()[1],fontsize=14,rotation='vertical',ha='center')

            if planes[plane]['index']==0:
                ax_masshist[row_index][planes[plane]['index']].set_ylabel(titles[row_index],fontsize=24)
        
        colourbar=fig_masshist.colorbar(imshows[-1],ax=ax_masshist[row_index],shrink=.95)
        colourbar.set_label('$log_{10}($Projected Density$)$ ($g/cm^2$)',fontsize=18)
        colourbar.ax.tick_params(labelsize=16)
        
        for ax in ax_masshist[row_index]:
            #ax.set_aspect('equal', adjustable='box')
            
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(18)
            #ax.yaxis.label.set_size(18)

            ax.set_xlim(-12,12)
            ax.set_ylim(-12,12)

            ax.tick_params(labelsize=16, axis='both',direction='in',top=True,right=True)
            ax.tick_params(axis='x', direction='in', pad=-20)
            ax.tick_params(axis='y', direction='in', pad=-30)

    plt.savefig(f'figures/{halo}/disrupteddensities.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def report_area_fracs():
    halos=['T1_Aug','T4_Aug','halo8']

    halo_names=['T1','T4','h8']

    halo_colours=['lightseagreen','deeppink','darkorange']

    bin_num=512

    fig_area_time,ax_area_time=plt.subplots(3,3,figsize=(10,8),sharey='row',sharex='col',constrained_layout=True)

    for halo in halos:
        planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}
        
        halo_index=halos.index(halo)

        for plane in planes:
        
            plane_index=planes[plane]['index']

            redshifts=np.load(f'halos/{halo}/area_fracs/{bin_num}px/redshifts.npy')
            area_frac_DLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/DLA.npy')[plane_index]  
            area_frac_subDLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/subDLA.npy')[plane_index]
            area_frac_LymanLimits=np.load(f'halos/{halo}/area_fracs/{bin_num}px/LymanLimit.npy')[plane_index]
            area_frac_lo_z_DLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/lo_z_DLA.npy')[plane_index]  
            area_frac_lo_z_subDLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/lo_z_subDLA.npy')[plane_index]
            area_frac_lo_z_LymanLimits=np.load(f'halos/{halo}/area_fracs/{bin_num}px/lo_z_LymanLimit.npy')[plane_index]

            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_DLAs),c='r',marker='d')
            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_subDLAs),c='b',marker='d')
            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_LymanLimits),c='g',marker='d')

            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_lo_z_DLAs),c='w',edgecolor='r',marker='d')
            ax_area_time[halo_index][plane_index].vlines(redshifts,np.array(area_frac_DLAs),np.array(area_frac_lo_z_DLAs),colors='r',ls='dashed',zorder=0,alpha=.5)

            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_lo_z_subDLAs),c='w',edgecolor='b',marker='d')
            ax_area_time[halo_index][plane_index].vlines(redshifts,np.array(area_frac_subDLAs),np.array(area_frac_lo_z_subDLAs),colors='b',ls='dashed',zorder=0,alpha=.5)

            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_lo_z_LymanLimits),c='w',edgecolor='g',marker='d')
            ax_area_time[halo_index][plane_index].vlines(redshifts,np.array(area_frac_LymanLimits),np.array(area_frac_lo_z_LymanLimits),colors='g',ls='dashed',zorder=0,alpha=.5)
            
            ax_area_time[halo_index][plane_index].invert_xaxis()
            ax_area_time[halo_index][plane_index].set_yscale('log')
            ax_area_time[halo_index][plane_index].set_facecolor(halo_colours[halo_index])
            ax_area_time[halo_index][plane_index].patch.set_alpha(0.1)

            ax_area_time[halo_index][plane_index].xaxis.label.set_size(16)
            ax_area_time[halo_index][plane_index].yaxis.label.set_size(16)
            ax_area_time[halo_index][plane_index].tick_params(labelsize=16)

            display_plane='{'+plane+'}'

            if halo_index==0:
                shift=0.27
            elif halo_index==1:
                shift=0.4
            else:
                shift=0.35

            ax_area_time[halo_index][plane_index].text(1.25,ax_area_time[halo_index][plane_index].get_ylim()[1]*shift,f'{halo_names[halo_index]}',fontsize=20,color=halo_colours[halo_index])
           

    
        ax_area_time[halo_index][1].tick_params(labelsize=16, axis='y',direction='in',which='both')
        ax_area_time[halo_index][2].tick_params(labelsize=16, axis='y',direction='in',which='both')

        ax_area_time[halo_index][1].set_ylabel('')
        ax_area_time[halo_index][2].set_ylabel('')

    ax_area_time[0][0].set_ylabel(r'Area Fraction',fontsize=16)
    ax_area_time[1][0].set_ylabel(r'Area Fraction',fontsize=16)
    ax_area_time[2][0].set_ylabel(r'Area Fraction',fontsize=16)

    ax_area_time[2][0].set_xlabel('Redshift ($z$)',fontsize=18)
    ax_area_time[2][1].set_xlabel('Redshift ($z$)',fontsize=18)
    ax_area_time[2][2].set_xlabel('Redshift ($z$)',fontsize=18)
    
    ax_area_time[0][0].set_title('$xy$ Plane',fontsize=18)
    ax_area_time[0][1].set_title('$xz$ Plane',fontsize=18)
    ax_area_time[0][2].set_title('$yz$ Plane',fontsize=18)

    ax_area_time[0][0].tick_params(labelsize=16, axis='x',direction='in',which='both')
    ax_area_time[1][0].tick_params(labelsize=16, axis='x',direction='in',which='both')

    ax_area_time[0][1].tick_params(labelsize=16, axis='x',direction='in',which='both') 
    ax_area_time[1][1].tick_params(labelsize=16, axis='x',direction='in',which='both')

    ax_area_time[0][2].tick_params(labelsize=16, axis='x',direction='in',which='both')
    ax_area_time[1][2].tick_params(labelsize=16, axis='x',direction='in',which='both')




    plt.savefig(f'figures/all_halos/area_fracs.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def report_px_num():

    halo='T1_Aug'
    snap_num=156
    plane='xy'


    bins=[128,256,512,1024]
    colours=['crimson','red','orangered','orange']
    panels=[(0,0),(1,0),(0,1),(1,1)]
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    fig_scatter, ax_scatter=plt.subplots(2,2,figsize=(10,8),constrained_layout=True,sharey='row',sharex='col')

    for bin_num in bins:

        col,row=panels[bins.index(bin_num)]

        nH_col=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/nH_col.npy')[planes[plane]['index']].flatten()
        mean_gz=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/mean_gz.npy')[planes[plane]['index']].flatten()
        masked_data={'DLA':{'param':'mean_gz'},'subDLA':{'param':'mean_gz'},'LymanLimit':{'param':'mean_gz'},'lo_z':{'param':'nH_col'}}
        
        for mask in masked_data:
            loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked/{mask}.npz')
            
            masked_data[mask]['data'] = np.ma.masked_array(loaded[f'{masked_data[mask]["param"]}_data'], mask=loaded[f'{masked_data[mask]["param"]}_mask'])

        lo_z_nH_col=masked_data['lo_z']['data'][planes[plane]['index']].compressed()
        DLA_mean_gz=masked_data['DLA']['data'][planes[plane]['index']].compressed()
        subDLA_mean_gz=masked_data['subDLA']['data'][planes[plane]['index']].compressed()
        LymanLimit_mean_gz=masked_data['LymanLimit']['data'][planes[plane]['index']].compressed()

        scatter=ax_scatter[row][col].scatter(nH_col,mean_gz,marker='x',c=colours[bins.index(bin_num)],s=2,alpha=.2,zorder=10)

        ax_scatter[row][col].text(2.5*10**9,6,r'$n_{px}=$'+str(bin_num),c=colours[bins.index(bin_num)],fontsize=18)

        xlims=[np.float64(0.8*10**9),np.float64(4*10**23)]
        ylims=[np.float64(0.3*10**-5),np.float64(2*10**1)]
        ax_scatter[row][col].set_yscale('log')
        ax_scatter[row][col].set_xscale('log')
        ax_scatter[row][col].set_xlim(xlims)
        ax_scatter[row][col].set_ylim(ylims)
        ax_scatter[row][col].tick_params(labelsize=16)

        ax_scatter[row][col].axvline(np.float64(10**20.3),c='r',ls='dashed')
        ax_scatter[row][col].fill_betweenx(np.array([10**-8,10**4]),np.float64(10**20.3),np.float64(10**24),color='r',alpha=.2)
        #ax_scatter[row][col].text(np.float64(0.2*10**22),np.float64(10**-2.5),'DLA',c='r',rotation=45,fontsize=14)

        ax_scatter[row][col].axvline(10**19,c='b',ls='dashed')
        ax_scatter[row][col].fill_betweenx(np.array([10**-8,10**4]),np.float64(10**19),np.float64(10**20.3),color='b',alpha=.2)
        #ax_scatter[row][col].text(np.float64(1.3*10**19),np.float64(10**-2.5),'Sub-DLA',c='b',rotation=45,fontsize=14)
        
        ax_scatter[row][col].axvline(10**17.2,c='g',ls='dashed')
        ax_scatter[row][col].fill_betweenx(np.array([10**-8,10**4]),np.float64(10**17.2),np.float64(10**19),color='g',alpha=.2)
        #ax_scatter[row][col].text(np.float64(6*10**17),np.float64(10**-2.5),'LLS',c='g',rotation=45,fontsize=14)

        ax_scatter[row][col].axhline(10**-3,color='blueviolet',ls='dashed',lw=2)
        #ax_scatter[row][col].text(np.float64(0.4*10**13),np.float64(1.2*10**-3),'Low Metallicity Threshold, $Z \leq 10^{-3}Z_{\odot}$',fontsize=18,color='blueviolet')

        ax_scatter[row][col].axhline(1,color='indigo',ls='dashed',lw=2)
        #ax_scatter[row][col].text(np.float64(0.8*10**16),np.float64(1.2),'Solar Metallicity, $Z \leq Z_{\odot}$',fontsize=18,color='indigo')

    ax_scatter[1][0].set_xlabel(f'Pixel Projected ${plane}$ Planar Number Density ($H_1^1/cm^2$)',fontsize=16)
    ax_scatter[1][0].xaxis.set_label_coords(1,-0.1)
    #ax_scatter[1][1].set_xlabel(f'Pixel Projected ${plane}$ Planar Number Density ($H_1^1/cm^2$)',fontsize=14)
    ax_scatter[1][0].set_ylabel('Solar-Relative Pixel-Mass-Weighted Mean Metallicity ($Z_\odot$)',fontsize=16)
    ax_scatter[1][0].yaxis.set_label_coords(-0.15,1)
    #ax_scatter[1][0].set_ylabel('Solar-Relative Pixel-Mass-Weighted Mean Metallicity ($Z_\odot$)',fontsize=12)

    ax_scatter[0][0].tick_params(axis='x', top=False, bottom=True, labeltop=False, labelbottom=False, direction='in') 
    ax_scatter[0][1].tick_params(axis='x', top=False, bottom=True, labeltop=False, labelbottom=False, direction='in') 

    ax_scatter[0][1].tick_params(axis='y', left=True, right=False, labelleft=False, labelright=False, direction='in',which='both') 
    ax_scatter[1][1].tick_params(axis='y', left=True, right=False, labelleft=False, labelright=False, direction='in',which='both')     
    
    plt.savefig('figures/banding.png',format="png",dpi=250,bbox_inches='tight')


    plt.show()

def T1_subfinds():
    if os.path.isdir(f'figures/T1_Aug/subfinds')!=True:
        os.makedirs(f'figures/T1_Aug/subfinds')

    fig_subfinds, ax_subfinds=plt.subplots(2,sharex=True)

    halo='T1_Aug'

    redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
    m200s=np.load(f'halos/{halo}/m200s/m200.npy')

    ax_subfinds[0].plot(redshifts,m200s,c='lightseagreen',zorder=3,lw=2)
    
    gas_zs=np.load(f'halos/{halo}/metallicity_behaviour/gas.npy')
    star_zs=np.load(f'halos/{halo}/metallicity_behaviour/stars.npy')

    if halo == 'T4_Aug':
        snap_0_gas_z=(0.0033585230966704242/z_sol+gas_zs[1])/2
        gas_zs[0]=snap_0_gas_z
        
    ax_subfinds[1].plot(redshifts,gas_zs,c='lightseagreen',zorder=3,lw=2)
    ax_subfinds[1].plot(redshifts,star_zs,c='lightseagreen',linestyle='dashdot',zorder=3,lw=2)

    for ax in ax_subfinds:
        ylims=ax.get_ylim()

        ax.fill_betweenx([ylims[0]/2,ylims[1]*2],5,2.208,facecolor='r',alpha=0.2)
        ax.axvline(2.208,color='w')
        ax.axvline(2.208,color='r',linestyle='dashed')

        ax.axvline(2.859,color='teal',linestyle='dotted')
        ax.axvline(2.494,color='teal',linestyle='dashed')

        ax.fill_betweenx([ylims[0]/2,ylims[1]*2],2.208,1.960,color=[1.0,0.65,0.65],facecolor='none',hatch='//')
        ax.axvline(1.960, color='darkviolet')

        ax.fill_betweenx([ylims[0]/2,ylims[1]*2],1.960,1.770,edgecolor='thistle',facecolor='none',hatch='//',alpha=.6)
        ax.axvline(1.770, color='w')
        ax.axvline(1.770, color='darkviolet',linestyle='dashed')

        ax.invert_xaxis()
        ax.set_xlim([4,1])
        ax.set_ylim(ylims)
        ax.set_yscale('log')

        ax.yaxis.set_tick_params(labelsize=14)
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())


    ax_subfinds[0].tick_params(axis='x', top=False, bottom=False, labeltop=False, labelbottom=False, direction='in') 
    ax_subfinds[1].tick_params(axis='x', bottom=True, labelbottom=True,direction='out',labelsize=14)


    
    ax_ticks = ax_subfinds[1].twiny()
    ax_ticks.set_xlim(ax_subfinds[1].get_xlim())

    ax_ticks.tick_params(axis='x', top=True, labeltop=False, direction='in')
    ax_ticks.spines['bottom'].set_visible(False)
    ax_ticks.spines['left'].set_visible(False)
    ax_ticks.spines['right'].set_visible(False)
    
    ax_subfinds[0].set_ylabel(r'$M_{200_{crit}} (M_\odot)$',fontsize=16)
    ax_subfinds[1].set_ylabel(r'Metallicity ($Z_\odot$)',fontsize=16)
    ax_subfinds[1].set_xlabel('Redshift (z)',fontsize=16)

    ax_subfinds[1].text(1.5,0.1,'$Z_{gas}$',c='lightseagreen',fontsize=16)
    ax_subfinds[1].text(1.3,0.165,'$Z_\star$',c='lightseagreen',fontsize=16)

    ax_subfinds[1].text(2.99,0.091,'Minor Merger', color='teal',rotation='vertical',fontsize=16,zorder=1)
    ax_subfinds[0].text(2.625,7.58*10**8,'Main Merger', color='teal',rotation='vertical',fontsize=16,zorder=1)

    ax_subfinds[1].text(2.14,0.13,'Merger-Induced SF',color='red',rotation='vertical',fontsize=16,zorder=1)

    ax_subfinds[1].text(1.92,0.155,'SNe Feedback',color='darkviolet',rotation='vertical',fontsize=16,zorder=1)

    plt.subplots_adjust(hspace=0)
    
    plt.savefig('figures/T1_Aug/subfinds/subfinds.pdf',format="pdf",dpi=250,bbox_inches='tight')
    plt.show()

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

    type_plot_info={'gas':{'cmap':'plasma','title':f'Gas'},'dm':{'cmap':'viridis','title':f'Dark Matter'},'stars':{'cmap':'magma','title':f'Stars'}}

    fig_masshist,ax_masshist=plt.subplots(len(req_types),3,figsize=(15,13),constrained_layout=True)
    
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
            ax_masshist[row_index][planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(8.5, 16),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))
            ax_masshist[row_index][planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(16, 8.5),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))                

            ax_masshist[row_index][planes[plane]['index']].text(12.25,9.5,planes[plane]['x_label'].split()[0],fontsize=18,ha='center')
            ax_masshist[row_index][planes[plane]['index']].text(9.5,12,planes[plane]['y_label'].split()[0],fontsize=18,rotation='vertical',ha='center')

            ax_masshist[row_index][planes[plane]['index']].text(11.75,6.5,planes[plane]['x_label'].split()[1],fontsize=18,ha='center')
            ax_masshist[row_index][planes[plane]['index']].text(7,10.25,planes[plane]['y_label'].split()[1],fontsize=18,rotation='vertical',ha='center')

            if planes[plane]['index']==0:
                ax_masshist[row_index][planes[plane]['index']].set_ylabel(type_plot_info[matter_type]['title'],fontsize=24)
        
        colourbar=fig_masshist.colorbar(imshows[-1],ax=ax_masshist[row_index],shrink=1)
        colourbar.set_label('$log_{10}($Projected Density$)$ ($g/cm^2$)',fontsize=18)
        colourbar.ax.tick_params(labelsize=16)
        
        for ax in ax_masshist[row_index]:
            #ax.set_aspect('equal', adjustable='box')
            
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(18)
            #ax.yaxis.label.set_size(18)

            ax.set_xlim(-37,37)
            ax.set_ylim(-37,37)

            ax.tick_params(labelsize=16, axis='both',direction='in',top=True,right=True)
            ax.tick_params(axis='x', direction='in', pad=-20)
            ax.tick_params(axis='y', direction='in', pad=-30)


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
        if matter_type=='dm':
            def logline(rho0,r,n):
                return(rho0*r**n)
            print(densities[0])
            #ax_radialdensity.plot(bin_centres[bin_centres<5],logline(3*10**-24,bin_centres[bin_centres<5],-1),ls='dashed',c='b',alpha=0.5)
            #ax_radialdensity.plot(bin_centres[bin_centres>0.5],logline(6*10**-24,bin_centres[bin_centres>0.5],-3),ls='dashed',c='b',alpha=0.5)
            #plt.text(2,logline(3*10**-24,4,-1)+0.1*10**-24,r'$\rho_{DM} \propto r^{-1}$',fontsize=12,c='b',rotation=-15)
            #plt.text(20,logline(6*10**-24,50,-3)+0.2*10**-28,r'$\rho_{DM} \propto r^{-3}$',fontsize=12,c='b',rotation=-38)
    
    ax_radialdensity.axvline(halo_r200,c='k',ls='dashed')
    plt.text(halo_r200+1, 10**-26, '$R_{200_{crit}}=$'+str(np.round(halo_r200,1))+'$ kpc$', fontsize=14,c='k',rotation='vertical')

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Radius ($kpc$)',fontsize=14)
    plt.ylabel('Spherical Radial Density ($g/cm^3$)',fontsize=14)
    plt.xlim(left=0.05)
    plt.ylim(top=10**-22)
    #plt.title(display_halo+', $z=$'+str(display_redshift),fontsize=16)
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
        req_dens=['nH_col']

    dens_plot_info={'total_mass':{'cmap':'plasma','title':f'{display_halo} Total Gas, z =${display_redshift}$, {bin_num} bins','cbar_label':'$log_{10}($Projected Density$)$ ($g/cm^2$)'},'hydrogen_mass':{'cmap':'viridis','title':f'{display_halo} Hydrogen Mass Fraction, z =${display_redshift}$, {bin_num} bins','cbar_label':'$log_{10}($Projected Density$)$ ($g/cm^2$)'},'nH_col':{'cmap':'cividis','title':f'{display_halo} Neutral Hydrogen Fraction, z =${display_redshift}$, {bin_num} bins','cbar_label':'$log_{10}($Column Density$)$ ($atoms/cm^2$)'}}

    fig_gashist,ax_gashist=plt.subplots(1,3,figsize=(15,5),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    for density in req_dens:
        loaded_data={f'{plane}':np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/{density}.npy')[planes[plane]['index']]for plane in planes}
        loaded_data['rel_pos']=read_raw_file(halo,'gas','rel_pos',snap_num=snap_num)
        
        vmin=np.log10(np.min([np.min(loaded_data[plane][loaded_data[plane]!=0.0]) for plane in planes])) #Find minimum (excluding zeros) and maximum projected density values across projection axes to normalise colour bars to
        vmax=np.log10(np.max([loaded_data[plane]for plane in planes]))

        extents=calc.get_extent(loaded_data['rel_pos'])
        plot_extents={plane:[extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc), extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)] for plane in planes}

        row_index=req_dens.index(density)

        imshows=[ax_gashist[planes[plane]['index']].imshow(np.log10(loaded_data[plane]),vmin=vmin,vmax=vmax,extent=plot_extents[plane],cmap=dens_plot_info[density]['cmap'],aspect='equal') for plane in planes]
        for plane in planes:
            ax_gashist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(8.5, 16),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))
            ax_gashist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(16, 8.5),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))                

            ax_gashist[planes[plane]['index']].text(12.25,9.5,planes[plane]['x_label'].split()[0],fontsize=18,ha='center')
            ax_gashist[planes[plane]['index']].text(9.5,12,planes[plane]['y_label'].split()[0],fontsize=18,rotation='vertical',ha='center')

            ax_gashist[planes[plane]['index']].text(11.75,6.5,planes[plane]['x_label'].split()[1],fontsize=18,ha='center')
            ax_gashist[planes[plane]['index']].text(7,10.25,planes[plane]['y_label'].split()[1],fontsize=18,rotation='vertical',ha='center')
        
        colourbar=fig_gashist.colorbar(imshows[-1],ax=ax_gashist,shrink=.9)

        colourbar.ax.axhline(20.3,c='r',lw=2)
        colourbar.ax.axhline(19,c='b',lw=2)
        colourbar.ax.axhline(17.2,c='g',lw=2)

        colourbar.ax.fill_betweenx([20.3,24],0,1,facecolor='none',edgecolor='red',hatch='///')
        colourbar.ax.fill_betweenx([19,20.3],0,1,facecolor='none',edgecolor='blue',hatch='///')
        colourbar.ax.fill_betweenx([17.2,19],0,1,facecolor='none',edgecolor='green',hatch='///')


        '''colourbar.ax.text(-2.2,20.3,'DLA',fontsize=16,c='r',rotation=45)
        colourbar.ax.text(-3.7,18,'Sub-DLA',fontsize=16,c='b',rotation=45)
        colourbar.ax.text(-2,17.5,'LLS',fontsize=16,c='g',rotation=45)'''
        colourbar.set_label(dens_plot_info[density]['cbar_label'],fontsize=18,ha='center')
        colourbar.ax.tick_params(labelsize=16)
        
        for ax in ax_gashist:
            #ax.set_aspect('equal', adjustable='box')
            
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(18)
            #ax.yaxis.label.set_size(18)

            ax.set_xlim(-37,37)
            ax.set_ylim(-37,37)

            ax.tick_params(labelsize=16, axis='both',direction='in',top=True,right=True)
            ax.tick_params(axis='x', direction='in', pad=-20)
            ax.tick_params(axis='y', direction='in', pad=-30)

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

    fig_gzhist,ax_gzhist=plt.subplots(1,3,figsize=(15,5),constrained_layout=True)
    
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
            ax_gzhist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(8.5, 16),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))
            ax_gzhist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(16, 8.5),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))                

            ax_gzhist[planes[plane]['index']].text(12.25,9.5,planes[plane]['x_label'].split()[0],fontsize=18,ha='center')
            ax_gzhist[planes[plane]['index']].text(9.5,12,planes[plane]['y_label'].split()[0],fontsize=18,rotation='vertical',ha='center')

            ax_gzhist[planes[plane]['index']].text(11.75,6.5,planes[plane]['x_label'].split()[1],fontsize=18,ha='center')
            ax_gzhist[planes[plane]['index']].text(7,10.25,planes[plane]['y_label'].split()[1],fontsize=18,rotation='vertical',ha='center')
        
    colourbar=fig_gzhist.colorbar(imshows[-1],ax=ax_gzhist,shrink=.9)
    colourbar.set_label('$log_{10}($Mass-Weighted Mean Solar-Relative Metallicity$)$ ($Z_\odot$)',fontsize=12)
    colourbar.ax.tick_params(labelsize=16)
        
    for ax in ax_gzhist:
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(18)
            #ax.yaxis.label.set_size(18)

            ax.set_xlim(-37,37)
            ax.set_ylim(-37,37)

            ax.tick_params(labelsize=16, axis='both',direction='in',top=True,right=True)
            ax.tick_params(axis='x', direction='in', pad=-20)
            ax.tick_params(axis='y', direction='in', pad=-30)

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

    plt.xlabel(f'Pixel Projected ${plane}$ Planar Number Density ($atoms/cm^2$)',fontsize=18)
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
    ax_scatter.text(np.float64(6*10**17),np.float64(10**-2.5),'LLS',c='g',rotation=45,fontsize=14)

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
        #plt.title(f'{display_halo} {plane} Projection, z =${display_redshift}$, {bin_num} bins',fontsize=20,pad=20)

        ax_projdenshist.tick_params(labelbottom=False,labelleft=False,labelright=True,labelsize=16)
        ax_projdenshist.yaxis.tick_right()
        ax_projdenshist.set_ylabel('$log_{10}(N_{px})$',fontsize=18)

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
        ax_gzhist.set_xlim([0,np.float64(10**5)])
        ax_gzhist.set_xscale('log')

        ax_gzhist.tick_params(labelleft=False,labelbottom=False,labeltop=True,labelsize=16)
        ax_gzhist.xaxis.tick_top()
        ax_gzhist.set_xlabel('$log_{10}(N_{px})$',fontsize=18)
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
    #plt.title(f'{display_halo}, {plane}-plane, {bin_num} bins',fontsize=20)
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

    fig_gashist,ax_gashist=plt.subplots(1,3,figsize=(15,5),constrained_layout=True)
    
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
            ax_gashist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(8.5, 16),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))
            ax_gashist[planes[plane]['index']].annotate('', xytext=(8.5, 8.5), xy=(16, 8.5),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))                

            ax_gashist[planes[plane]['index']].text(12.25,9.5,planes[plane]['x_label'].split()[0],fontsize=18,ha='center')
            ax_gashist[planes[plane]['index']].text(9.5,12,planes[plane]['y_label'].split()[0],fontsize=18,rotation='vertical',ha='center')

            ax_gashist[planes[plane]['index']].text(11.75,6.5,planes[plane]['x_label'].split()[1],fontsize=18,ha='center')
            ax_gashist[planes[plane]['index']].text(7,10.25,planes[plane]['y_label'].split()[1],fontsize=18,rotation='vertical',ha='center')
    
    colourbar=fig_gashist.colorbar(imshows[-1],ax=ax_gashist,shrink=.9)
    colourbar.set_label('$log_{10}($Projected Density$)$ ($g/cm^2$)',fontsize=18)
    colourbar.ax.tick_params(labelsize=16)
    
    for ax in ax_gashist:
            ax.set_box_aspect(1.0)
        
            ax.xaxis.label.set_size(18)
            #ax.yaxis.label.set_size(18)

            ax.set_xlim(-37,37)
            ax.set_ylim(-37,37)

            ax.tick_params(labelsize=16, axis='both',direction='in',top=True,right=True)
            ax.tick_params(axis='x', direction='in', pad=-20)
            ax.tick_params(axis='y', direction='in', pad=-30)

    if os.path.isdir(f'figures/{halo}/{bin_num}/contour_gas_densities')!=True:
        os.makedirs(f'figures/{halo}/{bin_num}/contour_gas_densities',exist_ok=True)

    plt.savefig(f'figures/{halo}/{bin_num}/contour_gas_densities/{snap_num}.pdf',format="pdf",dpi=250,bbox_inches='tight')
    
    plt.show()


def stellar_masses_redshift(halos):
    if os.path.isdir(f'figures/all_halos/stellar_masses')!=True:
        os.makedirs(f'figures/all_halos/stellar_masses')

    fig_masses,ax_masses=plt.subplots(figsize=(16,6))

    colours=['lightseagreen','deeppink','darkorange']
    markers=['D','D','D']

    nH_snaps={'T1_Aug':{'DLA':[[0,5],[6,7]],'subDLA':[[0,9],[11,13]],'LymanLimit':[[0,13]]},'halo8':{'DLA':[[0,1],[4,19],[26,29],[35,38]],'subDLA':[[0,38]],'LymanLimit':[[0,38]]},'T4_Aug':{'DLA':[[1,2],[3,4]],'subDLA':[[0,12]],'LymanLimit':[[0,13]]}}

    formation_events={'T1_Aug':{'Smaller Merger Begins':(2,'dashed'),'Small Merger Concludes, Main Merger Begins':(3,'dashed'),'Gas Merger Concludes':(5,'dashed'),'Overall Merger Concludes':(7,'dashed')},'T4_Aug':{'Merger Begins':(3,'dashed'),'Merger Concludes':(6,'dashed')},'halo8':{'Satellite In Halo':(15,'dotted'),'Central Disturbance':(19,'dotted'),'Disturbance Resolved':(21,'dotted'),'1st Merger Begins':(22,'dashed'),'1st Merger Concluding':(28,'dashed'),'2nd Merger In Progress':(31,'dashed'),'2nd Gas Merger Concludes':(32,'dashed'),'Tertiary Object In Halo':(33,'dotted'),'2nd Merger Concludes, 3rd Merger Begins':(34,'dashed'),'3rd Merger Concluding':(37,'dashed')}}


    for halo in halos:
        display_halo=halo.replace('_',' ')

        redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
        masses=np.load(f'halos/{halo}/stellar_masses/stellar_masses.npy')
          
        ax_masses.plot(redshifts,masses,c=colours[halos.index(halo)],zorder=1,lw=2)

        DLA_snaps=nH_snaps[halo]['DLA']
        subDLA_snaps=nH_snaps[halo]['subDLA']
        LLS_snaps=nH_snaps[halo]['LymanLimit']

        DLA_points=[ax_masses.scatter(redshifts[snaps[0]:snaps[1]],masses[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_masses.scatter(redshifts[snaps[0]:snaps[1]],masses[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_masses.scatter(redshifts[snaps[0]:snaps[1]],masses[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]

        events=[ax_masses.axvline(redshifts[formation_events[halo][event][0]],c=colours[halos.index(halo)],linestyle=formation_events[halo][event][1],alpha=0.5,zorder=0) for event in formation_events[halo]]

        event_labels=[ax_masses.text(redshifts[formation_events[halo][event][0]]+0.03,1.4*10**7,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor') if event=='Merger Begins' else ax_masses.text(redshifts[formation_events[halo][event][0]],1.3*10**7,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor')for event in formation_events[halo]]

        if halo=='halo8':
            shift=0.15
        else:
            shift=0.2
        plt.text(redshifts[0]+shift,masses[0],display_halo,fontsize=14,c=colours[halos.index(halo)])
    

    ax_masses.invert_xaxis()

    ax_masses.xaxis.label.set_size(18)
    ax_masses.yaxis.label.set_size(18)
    ax_masses.tick_params(labelsize=18)
    plt.xlabel('Redshift (z)',fontsize=18)
    plt.ylabel(r'Total Stellar Mass ($M_{\odot}$)',fontsize=18)
    plt.xlim([4.3,0.9])
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

    fig_zs,ax_zs=plt.subplots(figsize=(16,6))

    colours=['lightseagreen','deeppink','darkorange']
    markers=['D','D','D']

    nH_snaps={'T1_Aug':{'DLA':[[0,5],[6,7]],'subDLA':[[0,9],[11,13]],'LymanLimit':[[0,13]]},'halo8':{'DLA':[[0,1],[4,19],[26,29],[35,38]],'subDLA':[[0,38]],'LymanLimit':[[0,38]]},'T4_Aug':{'DLA':[[1,2],[3,4]],'subDLA':[[0,12]],'LymanLimit':[[0,13]]}}

    formation_events={'T1_Aug':{'Smaller Merger Begins':(2,'dashed'),'Small Merger Concludes, Main Merger Begins':(3,'dashed'),'Gas Merger Concludes':(5,'dashed'),'Overall Merger Concludes':(7,'dashed')},'T4_Aug':{'Merger Begins':(3,'dashed'),'Merger Concludes':(6,'dashed')},'halo8':{'Satellite In Halo':(15,'dotted'),'Central Disturbance':(19,'dotted'),'Disturbance Resolved':(21,'dotted'),'1st Merger Begins':(22,'dashed'),'1st Merger Concluding':(28,'dashed'),'2nd Merger In Progress':(31,'dashed'),'2nd Gas Merger Concludes':(32,'dashed'),'Tertiary Object In Halo':(33,'dotted'),'2nd Merger Concludes, 3rd Merger Begins':(34,'dashed'),'3rd Merger Concluding':(37,'dashed')}}

    for halo in halos:
        display_halo=halo.replace('_',' ')

        redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
        gas_zs=np.load(f'halos/{halo}/metallicity_behaviour/gas.npy')
        star_zs=np.load(f'halos/{halo}/metallicity_behaviour/stars.npy')

        if halo == 'T4_Aug':
            snap_0_gas_z=(0.0033585230966704242/z_sol+gas_zs[1])/2
            gas_zs[0]=snap_0_gas_z
          
        ax_zs.plot(redshifts,gas_zs,c=colours[halos.index(halo)],zorder=1,lw=2)
        ax_zs.plot(redshifts,star_zs,c=colours[halos.index(halo)],linestyle='dashdot',zorder=0,lw=2)

        print(gas_zs)

        DLA_snaps=nH_snaps[halo]['DLA']
        subDLA_snaps=nH_snaps[halo]['subDLA']
        LLS_snaps=nH_snaps[halo]['LymanLimit']

        DLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],gas_zs[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],gas_zs[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],gas_zs[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]

        DLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],star_zs[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],star_zs[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],star_zs[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]

        events=[ax_zs.axvline(redshifts[formation_events[halo][event][0]],c=colours[halos.index(halo)],linestyle=formation_events[halo][event][1],alpha=0.5,zorder=0) for event in formation_events[halo]]

        event_labels=[ax_zs.text(redshifts[formation_events[halo][event][0]]+0.03,0.245,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor') if event=='Merger Begins' else ax_zs.text(redshifts[formation_events[halo][event][0]],0.237,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor')for event in formation_events[halo]]

        if halo=='halo8':
            shift=0.31
        else:
            shift=0.35
        plt.text(redshifts[0]-0.03+shift,gas_zs[0],f'{display_halo} Gas',fontsize=14,c=colours[halos.index(halo)])
        plt.text(redshifts[0]+shift,star_zs[0],f'{display_halo} Stars',fontsize=14,c=colours[halos.index(halo)])

    ax_zs.invert_xaxis()

    ax_zs.xaxis.label.set_size(18)
    ax_zs.yaxis.label.set_size(18)
    ax_zs.tick_params(labelsize=18)
    plt.xlabel('Redshift (z)',fontsize=18)
    plt.ylabel(r'Metallicity ($Z_\odot$)',fontsize=18)
    plt.xlim([4.4,0.9])
    plt.yscale('log')

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

    fig_zs,ax_zs=plt.subplots(figsize=(16,6))

    colours=['lightseagreen','deeppink','darkorange']
    markers=['D','D','D']

    nH_snaps={'T1_Aug':{'DLA':[[0,5],[6,7]],'subDLA':[[0,9],[11,13]],'LymanLimit':[[0,13]]},'halo8':{'DLA':[[0,1],[4,19],[26,29],[35,38]],'subDLA':[[0,38]],'LymanLimit':[[0,38]]},'T4_Aug':{'DLA':[[1,2],[3,4]],'subDLA':[[0,12]],'LymanLimit':[[0,13]]}}

    formation_events={'T1_Aug':{'Smaller Merger Begins':(2,'dashed'),'Small Merger Concludes, Main Merger Begins':(3,'dashed'),'Gas Merger Concludes':(5,'dashed'),'Overall Merger Concludes':(7,'dashed')},'T4_Aug':{'Merger Begins':(3,'dashed'),'Merger Concludes':(6,'dashed')},'halo8':{'Satellite In Halo':(15,'dotted'),'Central Disturbance':(19,'dotted'),'Disturbance Resolved':(21,'dotted'),'1st Merger Begins':(22,'dashed'),'1st Merger Concluding':(28,'dashed'),'2nd Merger In Progress':(31,'dashed'),'2nd Gas Merger Concludes':(32,'dashed'),'Tertiary Object In Halo':(33,'dotted'),'2nd Merger Concludes, 3rd Merger Begins':(34,'dashed'),'3rd Merger Concluding':(37,'dashed')}}

    for halo in halos:
        display_halo=halo.replace('_',' ')

        redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
        m200s=np.load(f'halos/{halo}/m200s/m200.npy')

        DLA_snaps=nH_snaps[halo]['DLA']
        subDLA_snaps=nH_snaps[halo]['subDLA']
        LLS_snaps=nH_snaps[halo]['LymanLimit']

        DLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],m200s[snaps[0]:snaps[1]],c='r',marker=markers[halos.index(halo)],zorder=4,s=20) for snaps in DLA_snaps]
        subDLA_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],m200s[snaps[0]:snaps[1]],c='b',marker=markers[halos.index(halo)],zorder=3,s=20) for snaps in subDLA_snaps]
        LLS_points=[ax_zs.scatter(redshifts[snaps[0]:snaps[1]],m200s[snaps[0]:snaps[1]],c='g',marker=markers[halos.index(halo)],zorder=2,s=20) for snaps in LLS_snaps]

        events=[ax_zs.axvline(redshifts[formation_events[halo][event][0]],c=colours[halos.index(halo)],linestyle=formation_events[halo][event][1],alpha=0.5,zorder=0) for event in formation_events[halo]]

        event_labels=[ax_zs.text(redshifts[formation_events[halo][event][0]]+0.03,4.2*10**9,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor') if event=='Merger Begins' else ax_zs.text(redshifts[formation_events[halo][event][0]],4.1*10**9,event,c=colours[halos.index(halo)],rotation=45,rotation_mode='anchor')for event in formation_events[halo]]

        ax_zs.plot(redshifts,m200s,c=colours[halos.index(halo)],zorder=1,lw=2)
        if halo=='halo8':
            shift=0.15
        else:
            shift=0.2
        plt.text(redshifts[0]+shift,m200s[0],f'{display_halo}',fontsize=14,c=colours[halos.index(halo)])


    ax_zs.invert_xaxis()

    ax_zs.xaxis.label.set_size(18)
    ax_zs.yaxis.label.set_size(18)
    ax_zs.tick_params(labelsize=18)
    plt.xlabel('Redshift (z)',fontsize=18)
    plt.ylabel(r'$M_{200_{crit}} (M_\odot)$',fontsize=18)
    plt.xlim([4.3,0.9])
    plt.yscale('log')


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

    fig_weighted,ax_weighted=plt.subplots(figsize=(8,6))

    xlims=[np.float64(0.3*10**-3),np.float64(1)]

    plt.hist(unweighted_gz,edgecolor='blue',log=True,bins=np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),100),histtype='step')
    plt.hist(weighted_gz,color='dodgerblue',edgecolor='darkblue',log=True,bins=np.logspace(np.log10(xlims[0]),np.log10(xlims[1]),100),alpha=.4)
    
    plt.text(0.09,150,'Unweighted Mean',c='b',fontsize=16)
    

    plt.xlim(xlims)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Pixel Metallicity Value ($Z_\odot$)',fontsize=16)
    plt.ylabel('$log_{10}(N_{px})$',fontsize=16)

    ax_weighted.xaxis.label.set_size(18)
    ax_weighted.yaxis.label.set_size(18)
    ax_weighted.tick_params(labelsize=18)

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

    nh_masks=['DLA']

    fig_sightlineconts,ax_sightlineconts=plt.subplots(figsize=(7,7),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'}}

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
            
            raw_contours=[ax_sightlineconts.contour(rad_masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=raw_extents[plane],colors=colours[radii.index(radius)],zorder=3) for plane in planes if np.any(~rad_masks[col_dens]['data'][planes[plane]['index']].mask)]

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

            contours=[ax_sightlineconts.contour(rad_masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=extents[plane],colors=colours[radii.index(radius)],zorder=3,linewidths=5,alpha=0.7) for plane in planes if np.any(~rad_masks[col_dens]['data'][planes[plane]['index']].mask)]
            
        masks={'DLA':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.3],[1,1,1,0]]))},'subDLA':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.3],[1,1,1,0]]))},'LymanLimit':{'colour':'w','cmap':ListedColormap(np.array([[1,1,1,.3],[1,1,1,0]]))}}

        loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/512px/gas_only/masked/{col_dens}.npz') 
        masks[col_dens]['data'] = np.ma.masked_array(loaded[f'nH_col_data'], mask=loaded[f'nH_col_mask'])

        lo_z_loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/512px/gas_only/masked/{col_dens}.npz')
        masks[col_dens]['lo_z_data'] = np.ma.masked_array(lo_z_loaded[f'nH_col_data'], mask=lo_z_loaded[f'nH_col_mask'])

        raw_obj_contours=[ax_sightlineconts.contour(masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=obj_plot_extents[plane],colors=masks[col_dens]['colour']) for plane in planes if np.any(~masks[col_dens]['data'][planes[plane]['index']].mask)]

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

        adj_obj_contours=[ax_sightlineconts.contour(masks[col_dens]['data'][planes[plane]['index']][::-1, :].mask.astype(float),levels=[0.5],extent=shifted_obj_extents[plane],colors=masks[col_dens]['colour'],zorder=2,linewidths=5) for plane in planes if np.any(~masks[col_dens]['data'][planes[plane]['index']].mask)]
        adj_obj_fills=[ax_sightlineconts.imshow(masks[col_dens]['data'][planes[plane]['index']].mask.astype(float),extent=shifted_obj_extents[plane],cmap=masks[col_dens]['cmap'],zorder=1) for plane in planes if np.any(~masks[col_dens]['data'][planes[plane]['index']].mask)]  
        adj_obj_imshows=[ax_sightlineconts.imshow(np.log10(loaded_data[plane]),extent=shifted_obj_extents[plane],vmin=vmin,vmax=vmax,cmap=dens_plot_info['cmap'],aspect='equal',zorder=0) for plane in planes]
       
        for plane in planes:
            ax_sightlineconts.set_xlabel(planes[plane]['x_label'],fontsize=18)
            ax_sightlineconts.set_ylabel(planes[plane]['y_label'],fontsize=18)
            if planes[plane]['index']==1:
                ax_sightlineconts.set_title(col_dens,fontsize=20,pad=20)
    
    radii=[0.5,1,2,5,10]

    row_dims=[1,1.5,2.5]
    row_index=0

    ax_sightlineconts.set_xlim(-row_dims[row_index],row_dims[row_index])
    ax_sightlineconts.set_ylim(-row_dims[row_index],row_dims[row_index])
    ax_sightlineconts.set_aspect('equal', adjustable='box')
    ax_sightlineconts.set_aspect('equal', adjustable='box')
    ax_sightlineconts.set_box_aspect(1.0)
    ax_sightlineconts.xaxis.label.set_size(18)
    ax_sightlineconts.yaxis.label.set_size(18)

    ax_sightlineconts.yaxis.set_label_coords(-0.1, 0.5)

    ax_sightlineconts.tick_params(labelsize=16)

    ax_radscale=fig_sightlineconts.add_axes([0.965,0.117,0.12,0.833])
    


    ax_radscale.yaxis.tick_right()
    ax_radscale.set_ylim(0,5.5)
    ax_radscale.yaxis.set_label_position('right')
    ax_radscale.set_ylabel(r'Sample Annulus Radius ($R_{200_{crit}}$)',fontsize=18)
    ax_radscale.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax_radscale.xaxis.set_tick_params(labelsize=12)

    fractions=[0.5,1,2,5]
    labels=[r'$\frac{1}{2}R_{200_{crit}}$',r'$R_{200_{crit}}$',r'$2R_{200_{crit}}$',r'$5R_{200_{crit}}$']

    for frac in fractions:
        ax_radscale.axhline(frac,c=colours[fractions.index(frac)])
        ax_radscale.text(0.1,frac+0.1,labels[fractions.index(frac)],c=colours[fractions.index(frac)],fontsize=16)
    

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


def report_subfind_redshifts(halos):
    if os.path.isdir(f'figures/all_halos/subfinds')!=True:
        os.makedirs(f'figures/all_halos/subfinds')

    fig_subfinds, ax_subfinds=plt.subplots(1,2,figsize=(10,3),constrained_layout=True)

    colours=['lightseagreen','deeppink','darkorange']
    markers=['D','D','D']

    DLA_indexes={'T1_Aug':[0,5],'halo8':[8,19]}

    for halo in halos:

        redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
        m200s=np.load(f'halos/{halo}/m200s/m200.npy')

        ax_subfinds[0].plot(redshifts,m200s,c=colours[halos.index(halo)],zorder=1,lw=2)

        if halo!='T4_Aug':
            ax_subfinds[0].scatter(redshifts[DLA_indexes[halo][0]:DLA_indexes[halo][1]],m200s[DLA_indexes[halo][0]:DLA_indexes[halo][1]],c='r',marker='d')
        
        gas_zs=np.load(f'halos/{halo}/metallicity_behaviour/gas.npy')

        if halo == 'T4_Aug':
            ax_subfinds[1].plot(redshifts[0:2],[gas_zs[1],gas_zs[1]],c=colours[halos.index(halo)],linestyle='dashed',lw=2)
          
            ax_subfinds[1].plot(redshifts[1:],gas_zs[1:],c=colours[halos.index(halo)],zorder=1,lw=2)
        else:
            ax_subfinds[1].plot(redshifts,gas_zs,c=colours[halos.index(halo)],zorder=1,lw=2)
            ax_subfinds[1].scatter(redshifts[DLA_indexes[halo][0]:DLA_indexes[halo][1]],gas_zs[DLA_indexes[halo][0]:DLA_indexes[halo][1]],c='r',marker='d')


    for ax in ax_subfinds:
        ax.invert_xaxis()
        ax.set_xlim([4.1,1])
        ax.set_yscale('log')
        ax.tick_params(labelsize=16)
        ax.set_xlabel('Redshift (z)',fontsize=16)

    ax_subfinds[0].set_ylabel(r'$M_{200_{crit}} (M_\odot)$',fontsize=16)
    ax_subfinds[1].set_ylabel(r'Gas Metallicity ($Z_\odot$)',fontsize=16)

    ax_subfinds[0].text(2.1,2.75*10**9,'T1',c='lightseagreen',fontsize=16)
    ax_subfinds[0].text(1.5,9*10**8,'T4',c='deeppink',fontsize=16)
    ax_subfinds[0].text(1.7,1.5*10**9,'h8',c='darkorange',fontsize=16)

    ax_subfinds[1].text(1.6,0.08,'T1',c='lightseagreen',fontsize=16)
    ax_subfinds[1].text(2.7,0.15,'T4',c='deeppink',fontsize=16)
    ax_subfinds[1].text(3.2,0.04,'h8',c='darkorange',fontsize=16)

    
    plt.savefig('figures/all_halos/subfinds/subfinds.pdf',format="pdf",dpi=250,bbox_inches='tight')
    plt.show()

def report_T1_radial():
    
    snapshots=[152,155,156,157,164]
    colours=['red','teal','mediumblue','blueviolet','dodgerblue']

    fig_radialdensity,ax_radialdensity=plt.subplots()

    halo='T1_Aug'
    matter_type='gas'

    for snap_num in snapshots:
        loaded_data=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/radial_mass_density/{matter_type}.npy')
        
        bin_centres=loaded_data[0]
        densities=loaded_data[1]

        ax_radialdensity.plot(bin_centres,densities,c=colours[snapshots.index(snap_num)],alpha=0.7)
        #plt.text(bin_centres[0], densities[0], snap_num, fontsize=14,c=colours[snapshots.index(snap_num)])
       
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Radius ($kpc$)',fontsize=14)
    plt.ylabel('Spherical Radial Density ($g/cm^3$)',fontsize=14)
    plt.xlim(left=0.1,right=50)
    plt.ylim(top=10**-21)
    #plt.title(display_halo+', $z=$'+str(display_redshift),fontsize=16)
    ax_radialdensity.xaxis.set_tick_params(labelsize=12)
    ax_radialdensity.yaxis.set_tick_params(labelsize=12)

    #ax_radialdensity.text(32,0.15*10**-34,'T1',c='k',fontsize=20)

    pos=ax_radialdensity.get_position()

    ax_radialdensity.set_position([pos.x0, pos.y0, pos.width, pos.height*0.9])

    ax_redshifts=fig_radialdensity.add_axes([pos.x0, pos.y0+pos.height*0.9, pos.width, pos.height*0.1])

    redshifts=[4.008,2.494,2.208,1.960,0.997]
    labels=[r'$z\approx4$',r'Merger Begins',r'Gas Infall',r'SNe Feedback',r'$z\approx1$']
    shift=[-1.15,-3.1,-2.15,-3.05,-1.15]
    for redshift in redshifts:
        ax_redshifts.axvline(redshift,c=colours[redshifts.index(redshift)])
        ax_redshifts.text(redshift+0.07,shift[redshifts.index(redshift)],labels[redshifts.index(redshift)],c=colours[redshifts.index(redshift)],fontsize=14,rotation=300)
    
    ax_redshifts.invert_xaxis()
    ax_redshifts.xaxis.tick_top()
    ax_redshifts.xaxis.set_label_position('top')
    ax_redshifts.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    ax_redshifts.xaxis.set_tick_params(labelsize=12)
    ax_redshifts.set_xlabel('Redshift (z)',fontsize=14)
    ax_redshifts.set_xlim(4.2,0.75)


    plt.savefig('figures/T1_Aug/report_radials.pdf',format="pdf",dpi=250,bbox_inches='tight')


    plt.show()

def report_T4_densities():

    snap_num=152
    halo='T4_Aug'
    bin_num='512'
    matter_type='gas'
    title=f'$z=4.008$'
    plane='xz'


    type_plot_info={'gas':{'cmap':'plasma','title':f'Gas'},'dm':{'cmap':'viridis','title':f'Dark Matter'},'stars':{'cmap':'magma','title':f'Stars'}}

    fig_masshist,ax_masshist=plt.subplots(figsize=(6,6),constrained_layout=True)
    
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

   
    loaded_data={f'{plane}':np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/total_mass/{matter_type}.npy')[planes[plane]['index']]for plane in planes}
    loaded_data['rel_pos']=read_raw_file(halo,matter_type,'rel_pos',snap_num=snap_num)
        
    vmin=np.log10(np.min([np.min(loaded_data[plane][loaded_data[plane]!=0.0]) for plane in planes])) #Find minimum (excluding zeros) and maximum projected density values across projection axes to normalise colour bars to
    vmax=np.log10(np.max([loaded_data[plane]for plane in planes]))

    extents=calc.get_extent(loaded_data['rel_pos'])
    plot_extents={plane:[extents[planes[plane]['axes'][0]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][0]]['max'].to_value(units.kpc), extents[planes[plane]['axes'][1]]['min'].to_value(units.kpc),extents[planes[plane]['axes'][1]]['max'].to_value(units.kpc)] for plane in planes}

    imshows=ax_masshist.imshow(np.log10(loaded_data[plane]),extent=plot_extents[plane],vmin=vmin,vmax=vmax,cmap=type_plot_info[matter_type]['cmap'],aspect='equal')

    ax_masshist.annotate('', xytext=(5.5, 5.5), xy=(5.5, 11),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))
    ax_masshist.annotate('', xytext=(5.5, 5.5), xy=(11, 5.5),arrowprops=dict(arrowstyle="->",shrinkA=0, shrinkB=0))                

    ax_masshist.text(8.25,4.25,planes[plane]['x_label'].split()[0],fontsize=18,ha='center')
    ax_masshist.text(4.5,8,planes[plane]['y_label'].split()[0],fontsize=18,rotation='vertical',ha='center')

    ax_masshist.text(8.25,6.25,planes[plane]['x_label'].split()[1],fontsize=14,ha='center')
    ax_masshist.text(6.5,7.25,planes[plane]['y_label'].split()[1],fontsize=14,rotation='vertical',ha='center')


    ax_masshist.set_title(title,fontsize=24)
        
    colourbar=fig_masshist.colorbar(imshows,shrink=.75)
    colourbar.set_label('$log_{10}($Projected Density$)$ ($g/cm^2$)',fontsize=18)
    colourbar.ax.tick_params(labelsize=18)
            
    ax_masshist.set_box_aspect(1.0)

    ax_masshist.xaxis.label.set_size(18)
    #ax.yaxis.label.set_size(18)

    ax_masshist.set_xlim(-12,12)
    ax_masshist.set_ylim(-12,12)

    ax_masshist.tick_params(labelsize=18, axis='both',direction='in',top=True,right=True)
    ax_masshist.tick_params(axis='x', direction='in', pad=-20)
    ax_masshist.tick_params(axis='y', direction='in', pad=-30)

    plt.savefig(f'figures/{halo}/disrupteddensities.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def report_area_fracs():
    halos=['T1_Aug','T4_Aug','halo8']

    halo_names=['T1','T4','h8']

    halo_colours=['lightseagreen','deeppink','darkorange']

    bin_num=512

    fig_area_time,ax_area_time=plt.subplots(3,3,figsize=(10,8),sharey='row',sharex='col',constrained_layout=True)

    for halo in halos:
        planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}
        
        halo_index=halos.index(halo)

        for plane in planes:
        
            plane_index=planes[plane]['index']

            redshifts=np.load(f'halos/{halo}/area_fracs/{bin_num}px/redshifts.npy')
            area_frac_DLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/DLA.npy')[plane_index]  
            area_frac_subDLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/subDLA.npy')[plane_index]
            area_frac_LymanLimits=np.load(f'halos/{halo}/area_fracs/{bin_num}px/LymanLimit.npy')[plane_index]
            area_frac_lo_z_DLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/lo_z_DLA.npy')[plane_index]  
            area_frac_lo_z_subDLAs=np.load(f'halos/{halo}/area_fracs/{bin_num}px/lo_z_subDLA.npy')[plane_index]
            area_frac_lo_z_LymanLimits=np.load(f'halos/{halo}/area_fracs/{bin_num}px/lo_z_LymanLimit.npy')[plane_index]

            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_DLAs),c='r',marker='d')
            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_subDLAs),c='b',marker='d')
            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_LymanLimits),c='g',marker='d')

            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_lo_z_DLAs),c='w',edgecolor='r',marker='d')
            ax_area_time[halo_index][plane_index].vlines(redshifts,np.array(area_frac_DLAs),np.array(area_frac_lo_z_DLAs),colors='r',ls='dashed',zorder=0,alpha=.5)

            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_lo_z_subDLAs),c='w',edgecolor='b',marker='d')
            ax_area_time[halo_index][plane_index].vlines(redshifts,np.array(area_frac_subDLAs),np.array(area_frac_lo_z_subDLAs),colors='b',ls='dashed',zorder=0,alpha=.5)

            ax_area_time[halo_index][plane_index].scatter(redshifts,np.array(area_frac_lo_z_LymanLimits),c='w',edgecolor='g',marker='d')
            ax_area_time[halo_index][plane_index].vlines(redshifts,np.array(area_frac_LymanLimits),np.array(area_frac_lo_z_LymanLimits),colors='g',ls='dashed',zorder=0,alpha=.5)
            
            ax_area_time[halo_index][plane_index].invert_xaxis()
            ax_area_time[halo_index][plane_index].set_yscale('log')
            ax_area_time[halo_index][plane_index].set_facecolor(halo_colours[halo_index])
            ax_area_time[halo_index][plane_index].patch.set_alpha(0.1)

            ax_area_time[halo_index][plane_index].xaxis.label.set_size(16)
            ax_area_time[halo_index][plane_index].yaxis.label.set_size(16)
            ax_area_time[halo_index][plane_index].tick_params(labelsize=16)

            display_plane='{'+plane+'}'

            if halo_index==0:
                shift=0.27
            elif halo_index==1:
                shift=0.4
            else:
                shift=0.35

            ax_area_time[halo_index][plane_index].text(1.25,ax_area_time[halo_index][plane_index].get_ylim()[1]*shift,f'{halo_names[halo_index]}',fontsize=20,color=halo_colours[halo_index])
           

    
        ax_area_time[halo_index][1].tick_params(labelsize=16, axis='y',direction='in',which='both')
        ax_area_time[halo_index][2].tick_params(labelsize=16, axis='y',direction='in',which='both')

        ax_area_time[halo_index][1].set_ylabel('')
        ax_area_time[halo_index][2].set_ylabel('')

    ax_area_time[0][0].set_ylabel(r'Area Fraction',fontsize=16)
    ax_area_time[1][0].set_ylabel(r'Area Fraction',fontsize=16)
    ax_area_time[2][0].set_ylabel(r'Area Fraction',fontsize=16)

    ax_area_time[2][0].set_xlabel('Redshift ($z$)',fontsize=18)
    ax_area_time[2][1].set_xlabel('Redshift ($z$)',fontsize=18)
    ax_area_time[2][2].set_xlabel('Redshift ($z$)',fontsize=18)
    
    ax_area_time[0][0].set_title('$xy$ Plane',fontsize=18)
    ax_area_time[0][1].set_title('$xz$ Plane',fontsize=18)
    ax_area_time[0][2].set_title('$yz$ Plane',fontsize=18)

    ax_area_time[0][0].tick_params(labelsize=16, axis='x',direction='in',which='both')
    ax_area_time[1][0].tick_params(labelsize=16, axis='x',direction='in',which='both')

    ax_area_time[0][1].tick_params(labelsize=16, axis='x',direction='in',which='both') 
    ax_area_time[1][1].tick_params(labelsize=16, axis='x',direction='in',which='both')

    ax_area_time[0][2].tick_params(labelsize=16, axis='x',direction='in',which='both')
    ax_area_time[1][2].tick_params(labelsize=16, axis='x',direction='in',which='both')




    plt.savefig(f'figures/all_halos/area_fracs.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()

def report_px_num():

    halo='T1_Aug'
    snap_num=156
    plane='xy'


    bins=[128,256,512,1024]
    colours=['crimson','red','orangered','orange']
    panels=[(0,0),(1,0),(0,1),(1,1)]
    planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}

    fig_scatter, ax_scatter=plt.subplots(2,2,figsize=(10,8),constrained_layout=True,sharey='row',sharex='col')

    for bin_num in bins:

        col,row=panels[bins.index(bin_num)]

        nH_col=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/nH_col.npy')[planes[plane]['index']].flatten()
        mean_gz=np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/mean_gz.npy')[planes[plane]['index']].flatten()
        masked_data={'DLA':{'param':'mean_gz'},'subDLA':{'param':'mean_gz'},'LymanLimit':{'param':'mean_gz'},'lo_z':{'param':'nH_col'}}
        
        for mask in masked_data:
            loaded = np.load(f'/cosma/apps/durham/dc-coll7/halos/{halo}/{snap_num}/binned/{bin_num}px/gas_only/masked/{mask}.npz')
            
            masked_data[mask]['data'] = np.ma.masked_array(loaded[f'{masked_data[mask]["param"]}_data'], mask=loaded[f'{masked_data[mask]["param"]}_mask'])

        lo_z_nH_col=masked_data['lo_z']['data'][planes[plane]['index']].compressed()
        DLA_mean_gz=masked_data['DLA']['data'][planes[plane]['index']].compressed()
        subDLA_mean_gz=masked_data['subDLA']['data'][planes[plane]['index']].compressed()
        LymanLimit_mean_gz=masked_data['LymanLimit']['data'][planes[plane]['index']].compressed()

        scatter=ax_scatter[row][col].scatter(nH_col,mean_gz,marker='x',c=colours[bins.index(bin_num)],s=2,alpha=.2,zorder=10)

        ax_scatter[row][col].text(2.5*10**9,6,r'$n_{px}=$'+str(bin_num),c=colours[bins.index(bin_num)],fontsize=20)

        xlims=[np.float64(0.8*10**9),np.float64(4*10**23)]
        ylims=[np.float64(0.3*10**-5),np.float64(2*10**1)]
        ax_scatter[row][col].set_yscale('log')
        ax_scatter[row][col].set_xscale('log')
        ax_scatter[row][col].set_xlim(xlims)
        ax_scatter[row][col].set_ylim(ylims)
        ax_scatter[row][col].tick_params(labelsize=16)

        ax_scatter[row][col].axvline(np.float64(10**20.3),c='r',ls='dashed')
        ax_scatter[row][col].fill_betweenx(np.array([10**-8,10**4]),np.float64(10**20.3),np.float64(10**24),color='r',alpha=.2)
        #ax_scatter[row][col].text(np.float64(0.2*10**22),np.float64(10**-2.5),'DLA',c='r',rotation=45,fontsize=14)

        ax_scatter[row][col].axvline(10**19,c='b',ls='dashed')
        ax_scatter[row][col].fill_betweenx(np.array([10**-8,10**4]),np.float64(10**19),np.float64(10**20.3),color='b',alpha=.2)
        #ax_scatter[row][col].text(np.float64(1.3*10**19),np.float64(10**-2.5),'Sub-DLA',c='b',rotation=45,fontsize=14)
        
        ax_scatter[row][col].axvline(10**17.2,c='g',ls='dashed')
        ax_scatter[row][col].fill_betweenx(np.array([10**-8,10**4]),np.float64(10**17.2),np.float64(10**19),color='g',alpha=.2)
        #ax_scatter[row][col].text(np.float64(6*10**17),np.float64(10**-2.5),'LLS',c='g',rotation=45,fontsize=14)

        ax_scatter[row][col].axhline(10**-3,color='blueviolet',ls='dashed',lw=2)
        #ax_scatter[row][col].text(np.float64(0.4*10**13),np.float64(1.2*10**-3),'Low Metallicity Threshold, $Z \leq 10^{-3}Z_{\odot}$',fontsize=18,color='blueviolet')

        ax_scatter[row][col].axhline(1,color='indigo',ls='dashed',lw=2)
        #ax_scatter[row][col].text(np.float64(0.8*10**16),np.float64(1.2),'Solar Metallicity, $Z \leq Z_{\odot}$',fontsize=18,color='indigo')

    ax_scatter[1][0].set_xlabel(f'Pixel Projected ${plane}$ Planar Number Density ($H_1^1/cm^2$)',fontsize=18)
    ax_scatter[1][0].xaxis.set_label_coords(1,-0.1)
    #ax_scatter[1][1].set_xlabel(f'Pixel Projected ${plane}$ Planar Number Density ($H_1^1/cm^2$)',fontsize=14)
    ax_scatter[1][0].set_ylabel('Solar-Relative Pixel-Mass-Weighted Mean Metallicity ($Z_\odot$)',fontsize=18)
    ax_scatter[1][0].yaxis.set_label_coords(-0.15,1)
    #ax_scatter[1][0].set_ylabel('Solar-Relative Pixel-Mass-Weighted Mean Metallicity ($Z_\odot$)',fontsize=12)

    ax_scatter[0][0].tick_params(axis='x', top=False, bottom=True, labeltop=False, labelbottom=False, direction='in') 
    ax_scatter[0][1].tick_params(axis='x', top=False, bottom=True, labeltop=False, labelbottom=False, direction='in') 

    ax_scatter[0][1].tick_params(axis='y', left=True, right=False, labelleft=False, labelright=False, direction='in',which='both') 
    ax_scatter[1][1].tick_params(axis='y', left=True, right=False, labelleft=False, labelright=False, direction='in',which='both')     
    
    plt.savefig('figures/banding.png',format="png",dpi=250,bbox_inches='tight')


    plt.show()

def T1_subfinds():
    if os.path.isdir(f'figures/T1_Aug/subfinds')!=True:
        os.makedirs(f'figures/T1_Aug/subfinds')

    fig_subfinds, ax_subfinds=plt.subplots(2,sharex=True)

    halo='T1_Aug'

    redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
    m200s=np.load(f'halos/{halo}/m200s/m200.npy')

    ax_subfinds[0].plot(redshifts,m200s,c='lightseagreen',zorder=3,lw=2)
    
    gas_zs=np.load(f'halos/{halo}/metallicity_behaviour/gas.npy')
    star_zs=np.load(f'halos/{halo}/metallicity_behaviour/stars.npy')
        
    ax_subfinds[1].plot(redshifts,gas_zs,c='lightseagreen',zorder=3,lw=2)
    ax_subfinds[1].plot(redshifts,star_zs,c='lightseagreen',linestyle='dashdot',zorder=3,lw=2)

    for ax in ax_subfinds:
        ylims=ax.get_ylim()

        ax.fill_betweenx([ylims[0]/2,ylims[1]*2],5,2.208,facecolor='r',alpha=0.2)
        ax.axvline(2.208,color='w')
        ax.axvline(2.208,color='r',linestyle='dashed')

        ax.axvline(2.859,color='teal',linestyle='dotted')
        ax.axvline(2.494,color='teal',linestyle='dashed')

        ax.fill_betweenx([ylims[0]/2,ylims[1]*2],2.208,1.960,color=[1,0.7,0.6],facecolor='none',hatch='//')
        ax.axvline(1.960, color='darkviolet')

        ax.fill_betweenx([ylims[0]/2,ylims[1]*2],1.960,1.770,edgecolor='thistle',facecolor='none',hatch='//',alpha=.6)
        ax.axvline(1.770, color='w')
        ax.axvline(1.770, color='darkviolet',linestyle='dashed')

        ax.invert_xaxis()
        ax.set_xlim([4,1])
        ax.set_ylim(ylims)
        ax.set_yscale('log')

        ax.yaxis.set_tick_params(labelsize=14)
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())


    ax_subfinds[0].tick_params(axis='x', top=False, bottom=False, labeltop=False, labelbottom=False, direction='in') 
    ax_subfinds[1].tick_params(axis='x', bottom=True, labelbottom=True,direction='out',labelsize=14)
    
    ax_ticks = ax_subfinds[1].twiny()
    ax_ticks.set_xlim(ax_subfinds[1].get_xlim())

    ax_ticks.tick_params(axis='x', top=True, labeltop=False, direction='in')
    ax_ticks.spines['bottom'].set_visible(False)
    ax_ticks.spines['left'].set_visible(False)
    ax_ticks.spines['right'].set_visible(False)
    
    ax_subfinds[0].set_ylabel(r'$M_{200_{crit}} (M_\odot)$',fontsize=16)
    ax_subfinds[1].set_ylabel(r'Metallicity ($Z_\odot$)',fontsize=16)
    ax_subfinds[1].set_xlabel('Redshift (z)',fontsize=16)

    ax_subfinds[1].text(1.5,0.1,'$Z_{gas}$',c='lightseagreen',fontsize=16)
    ax_subfinds[1].text(1.3,0.165,'$Z_\star$',c='lightseagreen',fontsize=16)

    ax_subfinds[1].text(2.99,0.091,'Minor Merger', color='teal',rotation='vertical',fontsize=16,zorder=1)
    ax_subfinds[0].text(2.625,7.58*10**8,'Main Merger', color='teal',rotation='vertical',fontsize=16,zorder=1)

    ax_subfinds[1].text(2.14,0.13,'Merger-Induced SF',color='orangered',rotation='vertical',fontsize=16,zorder=1)

    ax_subfinds[1].text(1.92,0.155,'SNe Feedback',color='darkviolet',rotation='vertical',fontsize=16,zorder=1)

    plt.subplots_adjust(hspace=0)
    
    plt.savefig('figures/T1_Aug/subfinds/subfinds.pdf',format="pdf",dpi=250,bbox_inches='tight')
    plt.show()

def h8_subfinds():
    if os.path.isdir(f'figures/halo8/subfinds')!=True:
        os.makedirs(f'figures/halo8/subfinds')

    fig_subfinds, ax_subfinds=plt.subplots(2,sharex=True)

    halo='halo8'

    redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
    m200s=np.load(f'halos/{halo}/m200s/m200.npy')

    ax_subfinds[0].plot(redshifts,m200s,c='darkorange',zorder=10,lw=2)
    
    gas_zs=np.load(f'halos/{halo}/metallicity_behaviour/gas.npy')
    star_zs=np.load(f'halos/{halo}/metallicity_behaviour/stars.npy')
        
    ax_subfinds[1].plot(redshifts,gas_zs,c='darkorange',zorder=3,lw=2)
    ax_subfinds[1].plot(redshifts,star_zs,c='darkorange',linestyle='dashdot',zorder=3,lw=2)

    for ax in ax_subfinds:
        ylims=ax.get_ylim()

        ax.fill_betweenx([ylims[0]/2,ylims[1]*2],3.642,3.303,facecolor='red',alpha=0.2)
        ax.axvline(3.642,color='red')

        ax.fill_betweenx([ylims[0]/2,ylims[1]*2],3.642,3.427,edgecolor=[1,0.5,0.5],alpha=0.2,hatch='/',facecolor='none')
        ax.axvline(3.427,color='red',linestyle='dashed')

        ax.fill_betweenx([ylims[0]/2,ylims[1]*2],3.365,3.303,color=[1,0.7,0.6],facecolor='none',hatch='//',linewidth=0)
        ax.axvline(3.365,color='orangered',linestyle='dashed')
        ax.axvline(3.303, color='darkviolet')

        ax.fill_betweenx([ylims[0]/2,ylims[1]*2],3.303,3.203,edgecolor='thistle',facecolor='none',hatch='//',alpha=.6)
        ax.axvline(3.203, color='w')
        ax.axvline(3.203, color='darkviolet',linestyle='dashed')

        ax.axvline(3.163,color='orangered',linestyle='dotted')
        ax.axvline(3.104,color='orangered',linestyle='dotted')
        ax.axvline(2.544,color='orangered',linestyle='dotted')
        ax.axvline(1.496,color='orangered',linestyle='dotted')
        ax.axvline(1.185,color='orangered',linestyle='dotted')

        ax.invert_xaxis()
        ax.set_xlim([4,1])
        ax.set_ylim(ylims)
        ax.set_yscale('log')

        ax.yaxis.set_tick_params(labelsize=14)
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())


    ax_subfinds[0].tick_params(axis='x', top=False, bottom=False, labeltop=False, labelbottom=False, direction='in') 
    ax_subfinds[1].tick_params(axis='x', bottom=True, labelbottom=True,direction='out',labelsize=14)


    
    ax_ticks = ax_subfinds[1].twiny()
    ax_ticks.set_xlim(ax_subfinds[1].get_xlim())

    ax_ticks.tick_params(axis='x', top=True, labeltop=False, direction='in')
    ax_ticks.spines['bottom'].set_visible(False)
    ax_ticks.spines['left'].set_visible(False)
    ax_ticks.spines['right'].set_visible(False)
    
    ax_subfinds[0].set_ylabel(r'$M_{200_{crit}} (M_\odot)$',fontsize=16)
    ax_subfinds[1].set_ylabel(r'Metallicity ($Z_\odot$)',fontsize=16)
    ax_subfinds[1].set_xlabel('Redshift (z)',fontsize=16)

    ax_subfinds[1].text(1.45,0.02,'$Z_{gas}$',c='darkorange',fontsize=16)
    ax_subfinds[1].text(1.4,0.055,'$Z_\star$',c='darkorange',fontsize=16)

    ax_subfinds[1].text(3.6,0.04,'HI Accretion',rotation='vertical',color='red',fontsize=16,zorder=1)

    ax_subfinds[0].text(3.41,0.35*10**9,'Enriched SF',rotation='vertical',color='orangered',fontsize=16,zorder=3)
    ax_subfinds[0].fill_betweenx([0.33*10**9,2.52*10**9],3.42,3.33,facecolor=[1,0.8,0.8],zorder=2)

    ax_subfinds[0].text(3.25,0.33*10**9,'SNe Feedback',color='darkviolet',rotation='vertical',fontsize=16,zorder=3)
    ax_subfinds[0].fill_betweenx([0.3*10**9,3*10**9],3.24,3.15,facecolor='w',zorder=2)
    
    
    plt.subplots_adjust(hspace=0)
    
    plt.savefig('figures/halo8/subfinds/subfinds.pdf',format="pdf",dpi=250,bbox_inches='tight')
    plt.show()

def T4_subfinds():
    if os.path.isdir(f'figures/T4_Aug/subfinds')!=True:
        os.makedirs(f'figures/T4_Aug/subfinds')

    fig_subfinds, ax_subfinds=plt.subplots(2,sharex=True)

    halo='T4_Aug'

    redshifts=np.load(f'halos/{halo}/stellar_masses/redshifts.npy')
    m200s=np.load(f'halos/{halo}/m200s/m200.npy')

    ax_subfinds[0].plot(redshifts,m200s,c='deeppink',zorder=10,lw=2)
    
    gas_zs=np.load(f'halos/{halo}/metallicity_behaviour/gas.npy')
    star_zs=np.load(f'halos/{halo}/metallicity_behaviour/stars.npy')
        
    ax_subfinds[1].plot(redshifts[1:],gas_zs[1:],c='deeppink',zorder=3,lw=2)
    ax_subfinds[1].plot(redshifts,star_zs,c='deeppink',linestyle='dashdot',zorder=3,lw=2)
    ax_subfinds[1].plot(redshifts[0:2],[gas_zs[1],gas_zs[1]],zorder=2,lw=2,ls='dashed',c='deeppink')

    for ax in ax_subfinds:
        #ylims=ax.get_ylim()

        ax.invert_xaxis()
        ax.set_xlim([4,1])
        #ax.set_ylim(ylims)
        ax.set_yscale('log')

        ax.yaxis.set_tick_params(labelsize=14)
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())

        ax.axvline(2.494,color='mediumvioletred',linestyle='dotted')



    ax_subfinds[0].tick_params(axis='x', top=False, bottom=False, labeltop=False, labelbottom=False, direction='in') 
    ax_subfinds[1].tick_params(axis='x', bottom=True, labelbottom=True,direction='out',labelsize=14)


    
    ax_ticks = ax_subfinds[1].twiny()
    ax_ticks.set_xlim(ax_subfinds[1].get_xlim())

    ax_ticks.tick_params(axis='x', top=True, labeltop=False, direction='in')
    ax_ticks.spines['bottom'].set_visible(False)
    ax_ticks.spines['left'].set_visible(False)
    ax_ticks.spines['right'].set_visible(False)
    
    ax_subfinds[0].set_ylabel(r'$M_{200_{crit}} (M_\odot)$',fontsize=16)
    ax_subfinds[1].set_ylabel(r'Metallicity ($Z_\odot$)',fontsize=16)
    ax_subfinds[1].set_xlabel('Redshift (z)',fontsize=16)

    ax_subfinds[1].text(1.9,0.16,'$Z_{gas}$',c='deeppink',fontsize=16)
    ax_subfinds[1].text(1.4,0.137,'$Z_\star$',c='deeppink',fontsize=16)
    
    
    plt.subplots_adjust(hspace=0)
    
    plt.savefig('figures/T4_Aug/subfinds/subfinds.pdf',format="pdf",dpi=250,bbox_inches='tight')
    plt.show()

def report_mass_fracs():
    halos=['T1_Aug','halo8']

    halo_names=['T1','h8']

    halo_colours=['lightseagreen','darkorange']

    DLA_indexes={'T1_Aug':[0,5],'halo8':[8,19]}

    bin_num=512

    plane='xy'

    fig_area_time,ax_area_time=plt.subplots(figsize=(6,4),sharey='row',sharex='col',constrained_layout=True)


    for halo in halos:
        planes={'xy':{'index':0,'axes':['x','y'],'x_label':'$x$ ($kpc$)','y_label':'$y$ ($kpc$)'},'xz':{'index':1,'axes':['x','z'],'x_label':'$x$ ($kpc$)','y_label':'$z$ ($kpc$)'},'yz':{'index':2,'axes':['y','z'],'x_label':'$y$ ($kpc$)','y_label':'$z$ ($kpc$)'}}
        
        halo_index=halos.index(halo)
        
        plane_index=planes[plane]['index']

        redshifts=np.load(f'halos/{halo}/mass_fracs/{bin_num}px/redshifts.npy')
        area_frac_DLAs=np.load(f'halos/{halo}/mass_fracs/{bin_num}px/DLA.npy')[plane_index]  
        
        ax_area_time.plot(redshifts[DLA_indexes[halo][0]:DLA_indexes[halo][1]],np.array(area_frac_DLAs)[DLA_indexes[halo][0]:DLA_indexes[halo][1]],c=halo_colours[halo_index],alpha=0.5)
        ax_area_time.scatter(redshifts[DLA_indexes[halo][0]:DLA_indexes[halo][1]],np.array(area_frac_DLAs)[DLA_indexes[halo][0]:DLA_indexes[halo][1]],c=halo_colours[halo_index],marker='d',zorder=2)
       
        
    ax_area_time.invert_xaxis()
    ax_area_time.set_yscale('log')
    ax_area_time.set_xlim([3.65,3.29])
    ax_area_time.xaxis.label.set_size(16)
    ax_area_time.yaxis.label.set_size(16)
    ax_area_time.tick_params(labelsize=16)
    ax_area_time.tick_params(labelsize=16, axis='y',direction='in',which='both')
   

    ax_area_time.set_ylabel(r'$xy$ Mass Fraction',fontsize=18)
    ax_area_time.set_xlabel('Redshift ($z$)',fontsize=18)

    ax_area_time.text(3.6,0.04,'T1',c='lightseagreen',fontsize=18)
    ax_area_time.text(3.435,0.05,'h8',c='darkorange',fontsize=18)


    plt.savefig(f'figures/all_halos/mass_fracs.pdf',format="pdf",dpi=250,bbox_inches='tight')

    plt.show()