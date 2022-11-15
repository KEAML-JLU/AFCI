import random
import itertools
from src import config
from src.genetic import genetic
import numpy as np
import sys

#hopping
def hopping(population, parameters):
    J = parameters['jobs']
    jobNumber = len(J)
    jobsRange = range(1, jobNumber+1)
    sizeJobset1 = random.randint(0, jobNumber)

    jobset1 = random.sample(jobsRange, sizeJobset1)
    p1 = population[0]
    p2 = population[1]
    o1 = []
    p1kept = []
    for i in range(len(p1)):
        e = p1[i]
        if e in jobset1:
            o1.append(e)
        else:
            o1.append(-1)
            p1kept.append(e)

    o2 = []
    p2kept = []
    for i in range(len(p2)):
        e = p2[i]
        if e in jobset1:
            o2.append(e)
        else:
            o2.append(-1)
            p2kept.append(e)

    for i in range(len(o1)):
        if o1[i] == -1:
            o1[i] = p2kept.pop(0)

    for i in range(len(o2)):
        if o2[i] == -1:
            o2[i] = p1kept.pop(0)
    return population

#visiting
def visiting(population, parameters):
    newPop = []
    pop = population
    i = 0
    while i < len(population):
        (OS1, MS1) = population[i]
        (OS2, MS2) = population[i+1]

        if random.random() < config.pc:
            
            newPop.append((OS2, MS1))
            newPop.append((OS1, MS2))
        else:
            newPop.append((OS1, MS1))
            newPop.append((OS2, MS2))

        i = i + 2

    return pop

#releasing
def releasing(population, parameters):
    keptPopSize = int(config.pr * len(population))
    J = parameters['jobs']
    jobNumber = len(J)
    jobsRange = range(1, jobNumber+1)
    sizeJobset1 = random.randint(0, jobNumber)

    jobset1 = random.sample(jobsRange, sizeJobset1)
    p1 = population[0]
    p2 = population[1]
    o1 = []
    p1kept = []
    for i in range(len(p1)):
        e = p1[i]
        if e in jobset1:
            o1.append(e)
        else:
            o1.append(-1)
            p1kept.append(e)

    o2 = []
    p2kept = []
    for i in range(len(p2)):
        e = p2[i]
        if e in jobset1:
            o2.append(e)
        else:
            o2.append(-1)
            p2kept.append(e)

    for i in range(len(o1)):
        if o1[i] == -1:
            o1[i] = p2kept.pop(0)

    for i in range(len(o2)):
        if o2[i] == -1:
            o2[i] = p1kept.pop(0)
    return population

#assigning
def assigning(population, parameters):
    newPop = []
    i = 0
    while i < len(population):
        (OS1, MS1) = population[i]
        (OS2, MS2) = population[i+1]

        if random.random() < config.pc:
            
            newPop.append((OS2, MS1))
            newPop.append((OS1, MS2))
        else:
            newPop.append((OS1, MS1))
            newPop.append((OS2, MS2))

        i = i + 2
    return population

def choosing(start_ctr, duration, machine_jobs):
    max_duration_list = []
    max_duration = start_ctr + duration

    
    if machine_jobs:
        for job in machine_jobs:
            max_duration_list.append(job[3] + job[1])  # start + process time

        max_duration = max(max(max_duration_list), start_ctr) + duration

    machine_used = [True] * max_duration

    # Updating arr
    for job in machine_jobs:
        start = job[3]
        long = job[1]
        for k in range(start, start + long):
            machine_used[k] = False

    # Find the first place
    for k in range(start_ctr, len(machine_used)):
        if is_free(machine_used, k, duration):
            return k

def is_free(tab, start, duration):
    for k in range(start, start+duration):
        if not tab[k]:
            return False
    return True

def split_ms(pb_instance, ms):
    jobs = []
    current = 0
    for index, job in enumerate(pb_instance['jobs']):
        jobs.append(ms[current:current+len(job)])
        current += len(job)
    return jobs

def decode(pb_instance, os, ms):
    o = pb_instance['jobs']
    machine_operations = [[] for i in range(pb_instance['machinesNb'])]

    ms_s = split_ms(pb_instance, ms)  # machine for each operations

    indexes = [0] * len(ms_s)
    start_task_cstr = [0] * len(ms_s)

    # Iterating over OS to get task execution order and then checking in
    # MS to get the machine
    for job in os:
        index_machine = ms_s[job][indexes[job]]
        machine = o[job][indexes[job]][index_machine]['machine']
        prcTime = o[job][indexes[job]][index_machine]['processingTime']
        start_cstr = start_task_cstr[job]

        # Getting the first available place for the operation
        start = choosing(start_cstr, prcTime, machine_operations[machine - 1])
        name_task = "{}-{}".format(job, indexes[job]+1)

        machine_operations[machine - 1].append((name_task, prcTime, start_cstr, start))

        # Updating indexes (one for the current task for each job, one for the start constraint
        # for each job)
        indexes[job] += 1
        start_task_cstr[job] = (start + prcTime)

    return machine_operations