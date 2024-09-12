
# Created Jul 2024 by Cathleen Qiao

# This file is for reading radial data from a csv file. This was created 
# before I looked at the hfradarpy package. Thus, I was simply exploring
# what information the data could provide me. I recognized that there
# were timestamps and I could separate the data using them. I then
# plotted the longitutde and latitude for each time and tried other things
# such as removing the land flagged data. 



import pandas as pd
import matplotlib.pyplot as plt

# reading the csv file
data = pd.read_csv("./other-formatted-data/USF_MARA_May2024.csv", skiprows = [1], nrows = 50000)


# initial testing
''' 
data['longitude'] = pd.to_numeric(data['longitude'], errors = 'coerce')
data['latitude'] = pd.to_numeric(data['latitude'], errors = 'coerce')
data.dropna()

t1 = data[data['time'] == '2024-05-17T13:00:00Z']

t1.plot(kind='scatter', x='longitude', y='latitude')
plt.show()

print(data.describe())

data = data[data['VFLG'] != 128]

'''

# grouping the data by time for plotting
time = data.groupby('time')

#colors = data[data['VFLG'] != 128]

for l, t in time:
    t.plot(kind='scatter', x='longitude', y='latitude', c='VFLG', colormap='jet')
    plt.title(l)
    plt.show()
