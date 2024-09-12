
# Created Jul 2024 by Cathleen Qiao

# This section is for reading from a csv file. It provides basic information
# about the file such as the number of rows and the number of timestamps.


import pandas as pd

data = pd.read_csv("./other-formatted-data/USF_MARA_May2024.csv", skiprows=[1])

# data = pd.read_csv("USF_MARA_May2024.csv", skiprows=[1], dtype={'longitude': float})

unique_times = set(data['time'])

#print(unique_times)

print("Number of Rows: " + str(len(data)))
print("Numer of Times: " + str(len(unique_times)))

# data.info()


##########################################################################


# Created Jul 2024 by Cathleen Qiao

# This section is to provide more options when reading from a csv file. 
# You can print everything, print by line, and print separating the data
# into columns (for better visualization).


import csv

def printAll(r):
    for i, line in enumerate(r):
        print(f'line[{i}] = {line}')

def printLine(r):
    n = int(input("LINE: "))
    print(next(r))
    line = next((x for i, x  in enumerate(r) if i == n), None)
    print(line)

# print by line but separate the data by their columns
def printToCol(r):
    n = int(input("LINE: "))
    header = ''.join(next(r))
    line = ''.join(next((x for i, x  in enumerate(r) if i == n), None)) 
    for el, d in zip(header.split(','), line.split(',')):
        print(el + ": " + d)


with open("./other-formatted-data/USF_MARA_May2024.csv", "r") as f:
    reader = csv.reader(f, delimiter="\t")
    printToCol(reader)


