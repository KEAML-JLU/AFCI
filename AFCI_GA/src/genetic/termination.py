#!/usr/bin/env python
from src import config


def shouldTerminate(population, gen):
    return gen > config.maxGen
