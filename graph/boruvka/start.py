#!/usr/bin/python3
# Before start do it
# module add openmpi/4.0.0-icc slurm/15.08.1 java/jdk1.8.0
# HowTo sbatch -p test -n 5 --ntasks-per-node 1 --output=test.out run python start.py

import os
import subprocess
import time
import sys


def FileCheck(fn):
    try:
        f = open(fn, "r")
        f.close()
        return 1
    except IOError:
        return 0


if (len(sys.argv) != 5):
    print('Ti glupiy? Vvedi 4 parametra')
    print('"spark.sql.shuffle.partitions=?"')
    print('"--total-executor-cores", "?"')
    print('"--executor-cores", "?"')
    print('"--executor-memory", "?GB"')
    sys.exit(-1)

rank = os.environ["SLURM_PROCID"]
jobid = os.environ["SLURM_JOB_ID"]
n = int(os.environ["SLURM_JOB_NUM_NODES"])

if rank == "0":
    # Starting master
    pathToAddress = subprocess.check_output(["./spark/sbin/start-master.sh"]).split(" ")[4][:-1]
    time.sleep(10)

    # Writing master's address to file "jobid + /master_started"
    subprocess.check_output(["mkdir", jobid])
    file = open(pathToAddress, "r")
    s = file.read()
    file.close()
    masterAddress = s.split("\n")[14].split(" ")[8]

    file = open(jobid + "/master_started", "w")
    file.write(masterAddress)
    file.close()

    print("Master started! " + masterAddress + " " + pathToAddress)

    # Wait until all slaves started
    cnt = 1
    while cnt != n:
        cnt += FileCheck(jobid + "/slave_started_" + str(cnt))
        time.sleep(2)

    # Your spark submits
    # output = subprocess.check_output(["./spark/bin/spark-submit",
    #                                   "--class", "SparkPi",
    #                                   "--master", masterAddress,
    #                                   "--num-executors", "4",
    #                                   "--executor-cores", "6",
    #                                   "--executor-memory", "14GB",
    #                                   "sparkpi.jar", "24"])
    output = subprocess.check_output(["./spark/bin/spark-submit",
                                    "--conf", "spark.sql.shuffle.partitions=" + str(sys.argv[1]),
                                    "--class", "com.wrapper.BoruvkaAlgorithm",
                                    "--master", masterAddress,
                                    "--total-executor-cores", str(sys.argv[2]),
                                    "--executor-cores", str(sys.argv[3]),
                                    "--executor-memory", str(sys.argv[4]),
                                    "boruvka_2.11-0.1-SNAPSHOT.jar"])
    print(output)

    # Create finish file to signal slaves to stop
    subprocess.check_output(["touch", jobid + "/finish"])

    # Wait until all slaves stoped
    cnt = 1
    while cnt != n:
        cnt += FileCheck(jobid + "/slave_stoped_" + str(cnt))
        time.sleep(2)

    # Stoping master
    subprocess.check_output(["./spark/sbin/stop-master.sh"])
    time.sleep(10)
    subprocess.check_output(["touch", jobid + "/master_stoped"])
else:
    # Wait until master started and read master's address
    while(not FileCheck(jobid + "/master_started")):
        time.sleep(2)
    file = open(jobid + "/master_started", "r")
    masterAddress = file.read()
    file.close()

    # Starting slave
    subprocess.check_output(["./spark/sbin/start-slave.sh", masterAddress])
    time.sleep(10)

    # Create file to signal master that slave has started
    print("Slave started!")
    subprocess.check_output(["touch", jobid + "/slave_started_" + rank])

    # Wait for finish file
    while(not FileCheck(jobid + "/finish")):
        time.sleep(10)

    # Stoping slave and signal master
    subprocess.check_output(["./spark/sbin/stop-slave.sh"])
    time.sleep(10)
    subprocess.check_output(["touch", jobid + "/slave_stoped_" + rank])
