import sys
import random
import time
import math
from mpi4py import MPI as mpi
import numpy as np


def func(x, f):
    if f == 'sqr':
        return x * x
    elif f == 'sqrt':
        return math.sqrt(x)
    elif f == 'sin':
        return math.sin(x * math.pi)
    elif f == 'cube':
        return x * x * x
    elif f == 'hyp':
        return 1 / (1 + x)
    elif f == 'qcircle':
        return math.sqrt(1 - x * x)


def realInt(f):
    if f == 'sqr':
        return 1 / 3
    elif f == 'sqrt':
        return 2 / 3
    elif f == 'sin':
        return 2 / math.pi
    elif f == 'cube':
        return 1 / 4
    elif f == 'hyp':
        return math.log(2)
    elif f == 'qcircle':
        return math.pi / 4


if len(sys.argv) != 3:
    print('Check you comand string')
    sys.exit(-1)

comm = mpi.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
# print(size, rank)

n = 0.0
m = 0.0
random.seed(time.time() + rank)

# t1 = time.time()
t1 = mpi.Wtime()
for i in range(int(sys.argv[1])):
    if func(random.random(), sys.argv[2]) > random.random():
        m += 1
    n += 1

res = np.array(0.0)
r = np.array(m / n)
comm.Reduce(r, res, op=mpi.SUM, root=0)

# res = 0.0
# r = m / n
# res = comm.reduce(r, mpi.SUM, 0)

t2 = mpi.Wtime() - t1
# t2 = time.time() - t1

if rank == 0:
    # print(res)
    # print(r)
    print('Number of proc = ', size)
    print('Time =', t2)
    print('Real =', realInt(sys.argv[2]))
    print('Res = ', res / size)
