import h5py
import os
import pandas as pd
import csv



filename = "data_2023-02-22T17_29_37.489_Frohlich_T0.1_α3.0_N500_S20.csv"
datadir = "../data/"

dirt = os.path.join(datadir+datadir)

file = csv.reader(open(dirt))

for row in file:
    print(row)
    break

df = pd.read_csv(file, sep=';')