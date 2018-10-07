import sys
import random
import time
import math


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
        return math.ln(2)
    elif f == 'qcircle':
        return math.pi / 4


n = 0.0
m = 0.0

if len(sys.argv) != 3:
    print('Check you comand string')
    sys.exit(-1)

random.seed(time.time())

t1 = time.time()
for i in range(int(sys.argv[1])):
    if func(random.random(), sys.argv[2]) > random.random():
        m += 1
    n += 1

t2 = time.time() - t1

print('Number of proc = consistent')
print('Time =', t2)
print('Real =', realInt(sys.argv[2]))
print('Res = ', m / n)
