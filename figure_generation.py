from imports import *

def plot_projdens_xyz(xy,xz,yz,pos,ax,fig,matter_type,redshift,simulation): #Function to plot projected densities in 3 dimensions given xyz density data, position data and target MPL axes and figure

    types=['gas','dm','stars']
    maps=['plasma','viridis','magma']
    titles=[simulation+' Gas ($z=$'+str(redshift)+')',simulation+' Dark Matter ($z=$'+str(redshift)+')',simulation+' Stars ($z=$'+str(redshift)+')']
    colourbar_labels=['$log_{10}($Projected Number Density$)$ ($H_1^1/cm^2$)','$log_{10}($Projected Density$)$ ($g/cm^2$)','$log_{10}($Projected Density$)$ ($g/cm^2$)']

    colourmap=maps[types.index(matter_type)]
    title=titles[types.index(matter_type)]
    colourbar_label=colourbar_labels[types.index(matter_type)]

    vmin=np.log10(np.min([np.min(xy[xy!=0.0]),np.min(xz[xz!=0.0]),np.min(yz[yz!=0.0])])) #Find minimum (excluding zeros) and maximum projected density values across projection axes to normalise colour bars to
    vmax=np.log10(np.max([xy,xz,yz]))

    ax[0].imshow(np.log10(xy),extent=[np.min(pos[:,0]).value,np.max(pos[:,0]).value, np.min(pos[:,1]).value,np.max(pos[:,1]).value],vmin=vmin,vmax=vmax,aspect='equal',cmap=colourmap) #Plot imshows for first 2 dimensions
    ax[0].set_xlabel('$x$ ($kpc$)')
    ax[0].set_ylabel('$y$ ($kpc$)')
    
    ax[1].imshow(np.log10(xz),extent=[np.min(pos[:,0]).value,np.max(pos[:,0]).value, np.min(pos[:,2]).value,np.max(pos[:,2]).value],vmin=vmin,vmax=vmax,aspect='equal',cmap=colourmap)
    ax[1].set_xlabel('$x$ ($kpc$)')
    ax[1].set_ylabel('$z$ ($kpc$)')
    ax[1].set_title(title,fontsize=16,pad=20)
    
    colourbar=fig.colorbar(ax[2].imshow(np.log10(yz),extent=[np.min(pos[:,1]).value,np.max(pos[:,1]).value, np.min(pos[:,2]).value,np.max(pos[:,2]).value],vmin=vmin,vmax=vmax,aspect='equal',cmap=colourmap),ax=ax,pad=0.02,cmap=colourmap,shrink=0.75) #Plot final imshow and normalise colour bar for all subplots
    colourbar.set_label(colourbar_label,fontsize=12)
    ax[2].set_xlabel('$y$ ($kpc$)')
    ax[2].set_ylabel('$z$ ($kpc$)')
    
    for a in ax:
        a.set_aspect('equal', adjustable='box')
    