
# Created Aug 2024 by Cathleen Qiao

# This file is for getting the longitude and latitude of a point (P2) when 
# you are provided with the starting point (P1) longitude and latitude and 
# then the distance and bearing to the next point (P2). I wanted to find the 
# longitude and latitude of a point based on its range and bearing data (from the 
# antenna position) since it would help with my radial interpolation. However, 
# the accuracy was not high enough for me to use this method. 



# getting longitude and latitude from radial distance and bearing

from math import asin, atan2, degrees, radians, sin, cos

def get_point(lon1, lat1, r, b, R=6378.14):
    lat1 = radians(lat)
    lon1 = radians(lon)
    a = radians(b)
    lat2 = asin(sin(lat1) * cos(r / R) + cos(lat1) * sin(r / R) * cos(a))
    lon2 = lon1 + atan2(sin(a) * sin(r / R) * cos(lat1), cos(r / R) - sin(lat1) * sin(lat2))
    return (degrees(lon2), degrees(lat2))
    
lat = 24.7401333
lon = -80.9832833
r = 11.6497
b = 85.0

lon2, lat2 = get_point(lon, lat, r, b)

print(lon2, lat2)

# want -80.8695266
# want 24.7565425

#80.8685536
#24.7492557

# notes: only to 0.1 accuracy


