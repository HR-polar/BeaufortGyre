import numpy as npy
import matplotlib.pyplot as plt
import sys
import time
import os
import datetime
from netCDF4 import Dataset,num2date,date2num
import gyre_functions as gf

######################### IMPORTANT NOTES ######################
# This script assumes depth values are greater than zero. Be sure to multiply the depth array be -1 if this is the opposite!
# This script assumes -180 < longitude < 180. If it goes from 0-360, be sure to adjust for this!
# This script needs four input variables: SSH, latitude, longitude, and depth. All are required in netcdf format
# Parts that need to be edited or verified by the user can be found by #TOEDIT
######################### IMPORTANT NOTES ######################

######################### COMPUTATION SETTINGS #################
### Choose time to plot
s_year = YYYY #TOEDIT
e_year = YYYY #TOEDIT

timetype = 'monthly'; #TOEDIT ### Used in names
 
## This is the increment that will be iterated on for identifying closed contours
## smaller value = smaller increments, so takes longer, but less likely to miss small variations in contour edges so better for higher resolutions
increment_in = 0.05; incstr = '5cm' #TOEDIT

### locations to save outputs
model_name = '' #TOEDIT
outdir = ''     #TOEDIT
figdir = ''     #TOEDIT

plot_output = 1 # Set this to 1 for plotting a very basic map with BG contour. 0 if not

## Function to read in separate SSH files
def get_file_name_MODEL(year,month):

    ## If SSH fields are in different files, edit this function. If not, edit line further down in loop
    ## This function returns the SSH file for a given year and month
    ## It is used in the script when looping over files, so needs to be updated
    fname = fdir_monthly + FILENAME_YYYY_MM.nc #TOEDIT
    
    return fname

### set up locations of files to read - will vary depending on model
#TOEDIT BELOW
fdir_monthly = ''                           ## The location for the files with SSH 
thisvar = ''                                ## This is the variable name of SSH in the model
gdir = ''                                   ## The location for the file with bathymetry information
depthfile = ''                              ## The file containing bathymetry
depthvar = ''                               ## This is the variable name of bathymetry
latlondir = '';                             ## This is the location of latitude and longitude
latlonfile = ''                             ## This is the file containing latitude and longitude
latvar = ''; lonvar = ''                    ## Variable names for latitude and longitude
##########

lat = npy.squeeze(npy.array(Dataset(latlonfile).variables[latvar])) 
lon = npy.squeeze(npy.array(Dataset(latlonfile).variables[lonvar]))
depth = npy.squeeze(npy.array(Dataset(gdir+depthfile).variables[depthvar]))
land_mask = npy.zeros(npy.shape(lon))
land_mask[npy.abs(depth)>1e6] = 1 #TOEDIT ## This line should be modified so that it sets all land values to 1. Dependent on depth field/whether a land mask is already available
depth[land_mask==1] = 0 
######################### COMPUTATION SETTINGS DONE ############


##########################################################################
####Now loop to get gyre from different times

for c_year in range(s_year,e_year+1):
    
    for thism in range(0,num_months):
                
        incr = increment_in*1
 
        if timetype == 'monthly':
            timestr = 'y'+str(c_year)+'m'+str(thism+1).zfill(2) 
        elif timetype == 'yearly':
            timestr = 'y'+str(c_year)
        outname = outdir + 'BeaufortGyre_'+model_name+'_inc'+incstr+'_'+timetype+'_'+str(s_year)+'-'+str(e_year) ## TOEDIT
        print('------- OUTNAME IS '+outname)

        thisfile = get_file_name_MODEL(c_year,thism+1)
        
        print('about to start')

        if os.path.isfile(thisfile):
            
            print('starting')
            f_id = Dataset(thisfile)
            print('file '+thisfile)
            
            ssh_raw = npy.squeeze(npy.array(f_id.variables[thisvar])) ## if this has a time dimension, must modify this bit! It assumes a 2D variable
            ssh_raw[depth==0] = npy.nan
                
            start_time = time.time()
            [msk,lonmsk,latmsk,maxmsk]=gf.BG_compute(lon,lat,ssh_raw,depth,'SSH',incr,timestr,outname,[c_year,thism+1,15])
            print('saved '+outname+'...'+timestr)
            print(time.time() - start_time)

        if plot_output:
            ### plot to check quickly (no projection)
            plt.figure()
            px = plt.pcolormesh(ssh_raw,vmin=-0.2,vmax=0.2,cmap='RdYlBu_r')
            plt.colorbar(px)
            if npy.nansum(msk) > 0:
                msk_plot = msk*1
                msk_plot[npy.isnan(msk)] = 0
                plt.contour(msk_plot,colors='k')
            [r,c] = npy.nonzero(msk*ssh_raw==npy.nanmax(msk*ssh_raw))
            plt.plot(c,r,'k*')
            plt.contour(depth,[0,500,1000,1500,2000],colors='0.5',linewidths=0.5)
            plt.title(timestr)

            plt.savefig(figdir+'BG_from_'+model_name+'_'+timestr+'.png',bbox_inches='tight',dpi=100) #TOEDIT
            plt.close()

        ssh_raw = []
