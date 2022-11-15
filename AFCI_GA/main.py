#!/usr/bin/env python
import sys
import time
import pdb
import numpy as np
from src.utils import parser, gantt, AFCI
from src.genetic import encoding, decoding, genetic, termination
from src import config

import scipy.io as sio

# Beginning
if len(sys.argv) != 2:
    print("Usage: " + sys.argv[0] + " filename")
else:
    # Parameters Setting
    parameters = parser.parse(sys.argv[1])
    svg_name = sys.argv[1].split('\\')[-1].split('.')[0]
    # pdb.set_trace()

    

    # Initialize the pop
    pop = encoding.init_pop(parameters)
    gen,count,total_time = 1,0,0
    mm = [0]
    ave = []
    # pdb.set_trace()
    # Evaluate the pop
    # sP = np.zeros([0,480])
    while not termination.shouldTerminate(pop, gen):
        t0 = time.clock()
        # AFCI Process
        pop = AFCI.hopping(pop, parameters)
        pop = AFCI.visiting(pop, parameters)
        pop = AFCI.releasing(pop, parameters)
        pop = AFCI.assigning(pop, parameters)

        # Genetic Operators
        pop = genetic.selection(pop, parameters)
        pop = genetic.crossover(pop, parameters)
        pop = genetic.mutation (pop, parameters)
        t1 = time.clock()
        total_time += t1 - t0
        
        sortedPop = sorted(pop, key=lambda cpl: genetic.timeTaken(cpl, parameters))
        # if gen % 20 == 0 or gen == 1:
        #     savedsP = np.reshape(np.array(sortedPop),(500,-1))
        #     # pdb.set_trace()
        #     sP = np.append(sP, savedsP, axis=0)

        #savedsP = np.reshape(np.array(sortedPop),(500,-1))
        # pdb.set_trace()
        #sP = np.append(sP, savedsP, axis=0)

        # np.save('sPlocalbest.npy', sP, mode='a')
        # pdb.set_trace()
        m_p_m = genetic.timeTaken((sortedPop[0][0],sortedPop[0][1]),parameters)
        ave_p = genetic.timeTaken((sortedPop[-1][0],sortedPop[-1][1]),parameters)
        ave_ave = (ave_p + m_p_m) / 2
        if m_p_m == mm[-1]:
            count += 1
            
        else:
            count = 0
        mm.append(m_p_m)
        ave.append(ave_ave)
        gen = gen + 1
        if count == config.te:
            break
    # sP = np.reshape(sP, (-1,500,480))
    # datamat = 'F:\\Project\\PythonProject\\Python3\\flexible-job-shop-master\\sPGA.mat'
    # np.save('sPlocalbest.npy', sP)
    # matsP = sP.tolist()
    # sio.savemat(datamat, {'sP':matsP})
    sortedPop = sorted(pop, key=lambda cpl: genetic.timeTaken(cpl, parameters))
    
    # pdb.set_trace()
    
    
    print("Finished in {0:.2f}s".format(total_time))

    # Termination Criteria Satisfied ?
    gantt_data = decoding.translate_decoded_to_gantt(decoding.decode(parameters, sortedPop[0][0], sortedPop[0][1]))
    # m_p_m = genetic.timeTaken((sortedPop[0][0],sortedPop[0][1]),parameters)
    # pdb.set_trace()
    print(gen)
    print(mm)
    print(ave)
    if config.latex_export:
        gantt.export_latex(gantt_data)
    else:
        gantt.draw_chart(gantt_data,svg_name)
