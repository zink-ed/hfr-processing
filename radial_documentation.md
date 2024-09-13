# RADIAL_INTERPOLATION FILE – INTERPOLATION FUNCTION

This file contains the interpolation function for the radial data. Below is a detailed description of its functionality:

It takes in a Radial object (defined from the Radial file) and grabs the data from it (in the form of a dataframe). It then separates specific variables such as longitude, latitude, range cell, range, bearing, u velocity, v velocity, and radial velocity and converts them to numpy arrays.

- **1 = raw clean data (original)**

Then, for longitude, latitude, and bearing, dictionaries are created for each so that we can separate the data based on each range. So, the keys are the possible ranges and for each key / range, we have all the longitude, latitude, and bearing associated with it. This is simply to make it easier to interpolate based on nearest neighbors.

Then, going through each range, clusters are created based on bearing distance. The cluster starts with the first bearing. If the distance to the next bearing is more than 12 degrees or if the next bearing is the last bearing of the range, then the cluster ends. If the cluster has more than 2 points, then the longitude and latitude will be interpolated using a quadratic fit. The clusters are necessary so that the data is as accurate as possible, otherwise it will be interpolating long distances. Each interpolated value is added to the respective matrices and the interpolated matrix is updated. 

- **2 = interpolated longitude and latitude**

Now, to interpolate the velocities, we go through each cell of the interpolated matrix. If it has a flag = 2, then we can check to see if the value can be interpolated. We check points that are directly around it using the check_next function. It returns true if the cell has enough neighbor points that contain original data. It also appends the u and v velocities to a list where the interpolation function can then take the average of the points and update the cell for the matrices. If the cell doesn’t have enough direct neighbor points, then the get_grid function checks for farther neighbors and applies a bilinear interpolation formula to get the u and v velocities. In each case, the interpolated matrix is updated.

- **3 = interpolated u and v velocity based on neighbors**



