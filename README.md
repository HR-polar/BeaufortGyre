# BeaufortGyre
Scripts to analyse Beaufort Gyre in model simulations

read_and_compute_gyre_general.py reads in SSH, latitude, longitude, and depth variables and finds the largest closed contour 
of SSH (using gyre_functions.py) in the western Arctic Basin to identify the Beaufort Gyre. 

The code
1) reads in the monthly-averaged SSH field
2) finds the maximum SSH in the western Arctic (defined as a large box limited by longitude, latitude and also 1000m depth, 
to avoid capturing a maximum on the shelf). This is classed as the centre of the gyre 
3) an array is stored to hold a mask of the grid cells that fall within the gyre. The maximum SSH is the first grid cell to be added to this array
4) iterates out from this maximum, subtracting a given increment each time, checking new contours each time. 
For each iteration, it finds the list of points that comprise the edge of the new contour. I started with an initial increment of 1cm
- If this contour is closed, these new points within the contour are added to the gyre array and the algorithm continues for another 
iteration, subtracting the same increment again.
- If the contour is not closed, the algorithm tries a new increment which I set as 40% of the previous one, and checks this new contour. 
It keeps trying this increment adjustment until one is found that it can proceed with
5) Iterations continue until the size of the region within the contour does not change (or the increment gets sufficiently small).
To determine whether a contour is 'closed' or not, for each point in the given contour the algorithms check the surrounding points 
to check that they are not land. If they are, it defines it as 'not closed' and tries a smaller increment from the previous contour. 

Verification that the resulting contour was correct is done visually by looking at it on a map overlaying the SSH field.

User-specific variables to be changed are accompanied by #TOEDIT in read_and_compute_gyre_general.py. 
