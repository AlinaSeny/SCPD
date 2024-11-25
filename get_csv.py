import os
x = parallel_nums
y = data_size


d = {}
res_f = open('mpi_res.csv', 'w+')
files = os.listdir('./data/out/mpi')
for file in files:
    f = open('./data/out/mpi' + file, 'r')
    lines = f.readlines()
    for line in lines:
        if 'thread_num=' in line:
            tmp = findall(r'\d+', line)
            thread = tmp[0]
        if 'Nsize=' in line:
            tmp = findall(r'\d+', line)
            n = tmp[0]
        if 'avg_time=' in line:
            tmp = findall(r'\d+', line)
            t = tmp[0]
