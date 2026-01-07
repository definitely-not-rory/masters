import numpy as np 
import matplotlib.pyplot as plt 
import arepo_utils as ar
from scipy import stats as stats
from datetime import datetime
import nbimporter
import os
from matplotlib.patches import Circle
import astropy.units as units
import matplotlib as mpl
import matplotlib.font_manager as font_manager
import importlib

plt.rcParams["axes.formatter.use_mathtext"]=True
plt.rcParams['font.family']='serif'
cmfont = font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
plt.rcParams['font.serif']=cmfont.get_name()
plt.rcParams['mathtext.fontset']='cm'
plt.rcParams['axes.unicode_minus']=False

os.chdir('/cosma/home/durham/dc-coll7/masters') #Hard-resetting current working directory