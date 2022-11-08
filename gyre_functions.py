import numpy as npy
import matplotlib.pyplot as plt
import sys
import time
import os
import datetime
from netCDF4 import Dataset,num2date,date2num

def BG_compute(lon,lat,ssh_raw,depth,var_type,increment,timestr,outname,timearr,rm_landbarrier=0,do_npz=1):

  ## This function computes the largest closed contour in the Western Arctic basin
  ###### INPUTS:
  ##        lon:            longitude
  ##        lat:            latitude
  ##        ssh_raw:        the ssh field we're examining
  ##        depth:          the bathymetry (depth > 0) to check off-shelf regions 
  ##        var_type:       usually set to "SSH": identifies the variable type (has also been used for MSL in past)
  ##        increment:      the increment with which to iterate out from the maximum. Usually start at 10cm. 
  ##                        Higher resolution needs a smaller increment than lower resolution, because field varies more and so 
  ##                        larger increment may miss small features. But smaller increment takes longer
  ##        timestr:        a time string to use to name the .npz files (currently .nc does not have one, as all in one file)
  ##        outname:        the file name to output
  ##        rm_landbarrier  whether or not to use coastline as a valid edge of contour (e.g. set to 1 if using MSL, as atmospheric variable)
  ##        do_npz:         whether or not to save .npz files as well as .nc
  ###### OUTPUTS
  ##        mask_full:      the identified closed contour
  ##        lon:            longitude
  ##        lat:            latitude
  ##        BGcalcmin:      the ssh (or other variable) value at the edge of the closed contour

  ## Step 1: set up coastline to determine when the contour is no longer closed
  if rm_landbarrier==0: 
    land_arr = npy.nan*npy.ones([npy.shape(lon)[0],npy.shape(lon)[1]]); land_arr[npy.isnan(ssh_raw)] = 1;
    #First, reduce the array so that it takes less time
    lonmask = npy.nan*npy.ones(lon.shape);
    lonmask[lon<-80] = 1; lonmask[lon>140] = 1; lonmask[lat<68] = npy.nan;
    lonmask[(lon<20) & (lon > -130) & (lat < 70)] = npy.nan;        
    lonmask[(lon<20) & (lon > -120) & (lat < 75.5)] = npy.nan;
    lonmask[(lon<20) & (lon > -110) & (lat < 73)] = npy.nan;
    lonmask[(lon<20) & (lon > -100) & (lat < 80)] = npy.nan;
    lonmask[(lon<20) & (lon > -90) & (lat < 80) & (depth < 1000)] = npy.nan;
    ## extra
    lonmask[lat>80]  = npy.nan
    
    ssh_full = ssh_raw*lonmask;
    landmask = npy.zeros(lon.shape); landmask[ssh_raw==0] = 1; landmask[npy.isnan(ssh_raw)] = 1;

  else:
    land_arr = npy.zeros([npy.shape(lon)[0],npy.shape(lon)[1]]); land_arr[npy.isnan(ssh_raw)] = 1;
    lonmask = npy.ones(lon.shape)
    ssh_full = ssh_raw*lonmask

  ## Step 2: Identify off-shelf maximum    
  #This ensures that the maximum nonzero value off the shelf is found. Otherwise can get a high maxima near the coast
  ## Artificially force by depth field if rm_landbarrier isn't there 
  shelfmask = npy.ones(ssh_full.shape);
  if var_type == 'MSL':
    shelfmask[depth<0] = npy.nan; shelfmask[depth>=0] = 1;
  else:
    shelfmask[depth<3000] = npy.nan; shelfmask[depth>=3000] = 1; ## Here, make sure that depth array is in form of depths > 0!
  
  masked_shelf = ssh_full*shelfmask;
  masked_shelf[masked_shelf==0] = npy.nan;
  masked_shelf[npy.isnan(masked_shelf)] = npy.nan;
  maxarr = npy.nanmax(masked_shelf);
  print('max '+str(maxarr))
  maxarr_whole = maxarr.copy();

  if var_type == 'MSL':
    inc_min = increment/10000;
  else:
    inc_min = increment/1000; ## Was previously 100, for monthly values. Changed to 1000 for yearly. RERUN FOR MONTHLY
  
  #Here the array is reduced to a more manageable size
  ssh_xsum = npy.nansum(ssh_full*lonmask,axis=1);
  ssh_ysum = npy.nansum(ssh_full*lonmask,axis=0);
  xfind = npy.nonzero(ssh_xsum); x1 = npy.nanmin(xfind); x2 = npy.nanmax(xfind);
  yfind = npy.nonzero(ssh_ysum); y1 = npy.nanmin(yfind); y2 = npy.nanmax(yfind);
  
  ssh = ssh_full[x1:x2,y1:y2].copy(); dsmall = depth[x1:x2,y1:y2].copy();
  ssh_3000 = ssh.copy(); 
  if var_type != 'MSL':
    ssh_3000[dsmall<3000] = npy.nan

  mask_full = npy.zeros(ssh_full.shape);
  maskarr = npy.zeros(ssh.shape);
  land_arr = npy.zeros(ssh.shape); land_arr[npy.isnan(ssh)] = 1;
   
  if var_type=='FW': 
    land_arr[ssh==0] = 1; 
    land_arr[npy.isnan(ssh)] = 1; 

  ## Step 3: now loop over increments to find largest contour
  ## For each new increment, check that all cells in this new contour are a) not adjacent to land, and b) not higher than previous maximum
    
  ## LOOP #########################################
  #################################################
  if abs(maxarr) == 0:
      print('BG not found')
      all_met = 1;
  else:
      all_met = 0;
      for x in range (0,ssh.shape[0]):
          for y in range (0,ssh.shape[1]):
              if ssh_3000[x,y] >= maxarr_whole:
                  maskarr[x,y] = 1;

  size_of_old_mask = 0;
  size_of_new_mask = 0;
  while_loop = 0;
  cond = 1;
    
  ## We have the maximum value. Basically store checkarr coordinates and loop over it
  while all_met == 0:
      maskarr_new = maskarr.copy();
      reloop = 1;
      while_loop = while_loop + 1;

      ###################
      #Here the new maximum contour is found
      #Define a new edge array based on mask
      checkarr = fn_getEdge(maskarr_new);
      #Generate list of coordinates of new edge
      [cx,cy] = npy.where(checkarr==1);
      near_ocean = 1;

      length_of_carr = cx.shape[0];
      no_in_mask = npy.sum(maskarr_new+checkarr);
      inc_land_mask = maskarr_new+checkarr; inc_land_mask[land_arr==cond] = 0;
      if no_in_mask == npy.sum(inc_land_mask):
          looping = 1; end_of_loop = length_of_carr;
      else:
          looping = 1; end_of_loop = 1;
          near_ocean = 0; reloop = 0;

      while looping < end_of_loop:
          thisx = cx[looping];
          thisy = cy[looping];
          #First check that it's not land. If it is, exit the loop        
          if land_arr[thisx,thisy] != cond:
              #Need to check that this is next to the mask containing the maximum
              if ssh[thisx,thisy]>= maxarr:
                  maskarr_new[thisx,thisy] = 1;
                  #Now loop over surrounding cells
                  for yval in range (-1,2):
                      for xval in range (-1,2):
                          nexty = cy[looping]+yval; nextx = cx[looping]+xval;
                          #Check edges of domain
                          if nexty>=1: 
                              if nexty<ssh.shape[1]: 
                                  if nextx>=1: 
                                      if nextx<ssh.shape[0]:
                                          maxdim = npy.nanmax(ssh.shape);
                                          coords_1D = cy + 100*maxdim*maxdim*cx;
                                          maskarr_new[nextx,nexty] = 1;
                                          if nexty + 100*maxdim*maxdim*nextx not in coords_1D:
                                              cy=npy.append(cy,nexty);
                                              cx=npy.append(cx,nextx);
                                              end_of_loop = end_of_loop + 1;
          else:
              near_ocean = 0; looping = end_of_loop;

          looping = looping + 1;
          size_of_new_mask = npy.nansum(maskarr_new);
          maskarr_new_withland = maskarr_new.copy();
          maskarr_new_withland[land_arr==cond] = 0;
          size_of_new_mask_withland = npy.nansum(maskarr_new_withland);
      
      new_edges = fn_getEdge(maskarr_new)+maskarr_new;
      new_edges_withland = new_edges.copy(); new_edges_withland[land_arr==cond] = 0;
      #########################
      #check if its reached land or not
      if maxarr < npy.nanmin(ssh):
          all_met = 1; maskarr_out = maskarr.copy();
      
      if abs(increment) > abs(inc_min):        
          all_met = 0;
          if near_ocean == 0:
              maxarr = maxarr + increment;
              increment = increment/2;
              maxarr = maxarr - increment;
              maskarr_new = maskarr.copy();
              print('non-ocean cells. Try lower increment')
          elif npy.nansum(new_edges_withland) < npy.nansum(new_edges):
              maxarr = maxarr + increment;
              increment = increment/2;
              maxarr = maxarr - increment;
              maskarr_new = maskarr.copy();
              print('met a wall? reversing')
          elif size_of_new_mask > size_of_old_mask:
              print('continuing to increment out')
              maxarr = maxarr - increment;
              maskarr = maskarr_new.copy();
              size_of_old_mask = size_of_new_mask.copy();
          else:
              maxarr = maxarr + increment;
              increment = increment/2;
              maxarr = maxarr - increment;
              maskarr_new = maskarr.copy();
              print('other')
      else:
          all_met = 1; maskarr_out = maskarr.copy();
      
      print('next '+str(while_loop)+'max '+str(maxarr)+'no vals '+str(npy.sum(maskarr_new))+'inc '+str(increment))

  ############################
  mask_full[x1:x2,y1:y2] = maskarr_out.copy();

  BGcalcarr = mask_full*ssh_full;
  BGcalcarr[BGcalcarr==0] = npy.nan
  BGcalcmin = npy.nanmin(BGcalcarr)
  print('saving in '+outname)
  thisname_notime = outname+'fromSSH'
  if rm_landbarrier:
    thisname = thisname_notime+timestr+'nolandbarrier'
  else:
    thisname = thisname_notime+timestr
    
  if do_npz:  
    npy.savez(thisname+'.npz',BGlon = lon, BGlat = lat, BGmask = mask_full,BGmin = BGcalcmin)

  ## netcdf file output
  mask_nan = mask_full*1; mask_full[mask_full!=1] = npy.nan
  ### set up netcdf to save into
  outnetcdf = thisname_notime+'.nc'
  timevalue = datetime.datetime(timearr[0],timearr[1],timearr[2])
  time_unit_out = "seconds since 1979-01-01 00:00:00"

  if os.path.isfile(outnetcdf)==0:
      print('writing')
      ncfile = Dataset(outnetcdf,mode='w',format='NETCDF4_CLASSIC') 

      y = ncfile.createDimension('y',npy.shape(lon)[0])
      x = ncfile.createDimension('x',npy.shape(lon)[1])
      time = ncfile.createDimension('time',None)

      lat_var = ncfile.createVariable('BGlat', npy.float32, ('y','x'))
      lat_var.units = 'degrees_north'
      lat_var.long_name = 'latitude'
      lon_var = ncfile.createVariable('BGlon', npy.float32, ('y','x'))
      lon_var.units = 'degrees_east'
      lon_var.long_name = 'longitude'
      time = ncfile.createVariable('time', npy.float64, ('time',))
      time.units = time_unit_out
      time.long_name = 'time'

      lat_var[:,:] = lat
      lon_var[:,:] = lon
      time[:] = date2num(timevalue,time_unit_out)  

      msk_var = ncfile.createVariable('BGmask',npy.float64,('time','y','x')) 
      msk_var.units = '[]' 
      msk_var.standard_name = 'Beaufort_Gyre_mask' 
      BG_max = ncfile.createVariable('BGmax',npy.float64,('time')) 
      BG_max.units = 'metres' 
      BG_max.standard_name = 'Maximum_ssh_in_gyre'
      BG_min = ncfile.createVariable('BGmin',npy.float64,('time')) 
      BG_min.units = 'metres' 
      BG_min.standard_name = 'Minimum_ssh_in_gyre'

      msk_var[0,:,:] = mask_full
      BG_max[0] = npy.nanmax(mask_full*ssh_full)
      BG_min[0] = BGcalcmin

      ncfile.close()
  else:
      ## append to existing netcdf
      f = Dataset(outnetcdf,'a')
      ## first, append the time
      appendvar = f.variables['time']; 
      len_file = len(appendvar)
      if date2num(timevalue,time_unit_out) in appendvar:
            print('there is already an entry for this time - not appending')
        
      else:

          print('appending')
      
          appendvar[len_file] = date2num(timevalue,time_unit_out)

          ## now, save each variable
          appendvar = f.variables['BGmask']; appendvar[len_file,:,:] = mask_full
          appendvar = f.variables['BGmin']; appendvar[len_file] = BGcalcmin
          appendvar = f.variables['BGmax']; appendvar[len_file] = npy.nanmax(mask_full*ssh_full)

      f.close()

  return mask_full, lon, lat, BGcalcmin ### returns the full mask, longitude, latitude, and SSH at gyre edge

def fn_getEdge(oldarr):

    ## This function identifies the edge of a contour by looking at the four adjacent cells
    
    #Oldarr is a mask. Newarr finds coordinates next to ones
    newarr = npy.zeros(oldarr.shape);    
    
    dx = npy.diff(oldarr,axis=0);
    dy = npy.diff(oldarr,axis=1);
    y_offset_top = npy.zeros(dy.shape);
    y_offset_bottom = npy.zeros(dy.shape);
    x_offset_left = npy.zeros(dx.shape);
    x_offset_right = npy.zeros(dx.shape);
        
    y_offset_top[dy==1] = 1; 
    y_offset_bottom[dy==-1] = 1;
    x_offset_left[dx==1] = 1;
    x_offset_right[dx==-1] = 1;
        
    ## putting it into new mask
    newmask_y = npy.zeros(newarr.shape);
    newmask_x = npy.zeros(newarr.shape);
    newmask_y[:,0] = y_offset_top[:,0];
    newmask_y[:,-1] = y_offset_bottom[:,-1];
    newmask_y[:,1:-2] = y_offset_top[:,1:-1] + y_offset_bottom[:,0:-2];
    newmask_x[0,:] = x_offset_left[0,:];
    newmask_x[-1,:] = x_offset_right[-1,:];
    newmask_x[1:-2,:] = x_offset_left[1:-1,:] + x_offset_right[0:-2,:];
       
    newarr[newmask_x+newmask_y > 0] = 1; 
        
    return newarr
