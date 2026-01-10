import numpy as np #Modules And Packages Required
import matplotlib.pyplot as plt 
import arepo_utils as ar
from scipy import stats as stats
import os
from matplotlib.patches import Circle
import astropy.units as units
import matplotlib as mpl
from matplotlib import cm
import matplotlib.font_manager as font_manager
import importlib
import sys
from matplotlib.colors import ListedColormap

plt.rcParams["axes.formatter.use_mathtext"]=True #Hard-codes figure labels, fonts etc.
plt.rcParams['font.family']='serif'
cmfont = font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
plt.rcParams['font.serif']=cmfont.get_name()
plt.rcParams['mathtext.fontset']='cm'
plt.rcParams['axes.unicode_minus']=False

os.chdir('/cosma/home/durham/dc-coll7/masters') #Hard-resetting current working directory

z_sol=0.0196
mass_H_atom=1.6735575*10**-24