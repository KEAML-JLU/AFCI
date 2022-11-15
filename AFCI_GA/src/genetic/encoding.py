#!/usr/bin/env python

# This module creates a population of random OS and MS chromosomes

import random
from src import config


def gen_OS(parameters):
    jobs = parameters['jobs']

    OS = []
    i = 0
    for job in jobs:
        for op in job:
            OS.append(i)
        i = i+1

    random.shuffle(OS)

    return OS


def gen_MS(parameters):
    jobs = parameters['jobs']

    MS = []
    for job in jobs:
        for op in job:
            randomMachine = random.randint(0, len(op)-1)
            MS.append(randomMachine)

    return MS


def init_pop(parameters):
    gen1 = []

    for i in range(config.popSize):
        OS = gen_OS(parameters)
        MS = gen_MS(parameters)
        gen1.append((OS, MS))
        # genl = init_group(gen1, parameters)

    return gen1

# def init_group(Population, parameters):


#     return genl