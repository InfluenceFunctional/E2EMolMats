import numpy as np

data_path = r'D:\crystal_datasets\phonons\dynmat.dat'

lines = []
with open(data_path, 'r') as f:
    for line in f:
        lines.append(f.readline().split(' '))

lines = np.array(lines).astype('float')



