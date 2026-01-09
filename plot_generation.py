from imports import *
from halo_readers import get_redshift, get_snap_num

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

def proj_mass_density(halo,**kwargs):
    if 'snap_num' not in kwargs:
        if 'redshift' in kwargs:
            target_redshift=kwargs['redshift']
            snap_num,snap_redshift=get_snap_num(halo,target_redshift)
        else:
            sys.exit('Please provide either a target redshift (\"redshift=X\") or snapshot number (\"snap_num=XXX\")')
    else:
        snap_num=kwargs['snap_num']
        snap_redshift=get_redshift(halo,snap_num)
