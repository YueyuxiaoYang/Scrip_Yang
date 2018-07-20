# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 22:40:48 2018

@author: kevin1024
"""
#    This file is part of DEAP.
#
#    DEAP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    DEAP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with DEAP. If not, see <http://www.gnu.org/licenses/>.

import random

import numpy 

from deap import algorithms
from deap import base
from deap import creator
from deap import tools
from getSignal import *


def evalOneMax(individual,sum_signal_art,Parameters):
    return sum(pow((sumSignal(np.where(individual==1)[0],Parameters)-sum_signal_art),2)),

def cxTwoPointCopy(ind1, ind2):
    """Execute a two points crossover with copy on the input individuals. The
    copy is required because the slicing in numpy returns a view of the data,
    which leads to a self overwritting in the swap operation. It prevents
    ::
    
        >>> import numpy
        >>> a = numpy.array((1,2,3,4))
        >>> b = numpy.array((5.6.7.8))
        >>> a[1:3], b[1:3] = b[1:3], a[1:3]
        >>> print(a)
        [1 6 7 4]
        >>> print(b)
        [5 6 7 8]
    """
    size = len(ind1)
    cxpoint1 = random.randint(1, size)
    cxpoint2 = random.randint(1, size - 1)
    if cxpoint2 >= cxpoint1:
        cxpoint2 += 1
    else: # Swap the two cx points
        cxpoint1, cxpoint2 = cxpoint2, cxpoint1

    ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] \
        = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()
        
    return ind1, ind2

def generate_individual(IND_SIZE,Nbr_poly_estimate=150):
    Poly_position = np.random.choice(range(int(IND_SIZE)),size=Nbr_poly_estimate,replace=False)
    Pattern_poly = np.zeros(int(IND_SIZE))
    Pattern_poly[Poly_position] = 1
    return Pattern_poly

    
'''
Genetic Algorithm 
input:
    sum_signal_art: objective function which you want to be fitted, here is the
                    sum signals of artificial/experiment data
    Parameters
    IND_SIZE: length of polymerase position vector(Pattern_poly)
    population_Nbr: number of population in GA
    
'''
def GA(sum_signal_art,Parameters,IND_SIZE = 4,Population_Nbr=50,Max_gen=20,Nbr_poly_estimate=150):
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", numpy.ndarray, fitness=creator.FitnessMin)
    
    toolbox = base.Toolbox()
    
    #toolbox.register("attr_bool", random.randint, 0, 1)
    toolbox.register("generate_individual", generate_individual, IND_SIZE,Nbr_poly_estimate)
    toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.generate_individual)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)                                                                                      
    
    toolbox.register("evaluate", evalOneMax,sum_signal_art=sum_signal_art,Parameters=Parameters)
    toolbox.register("mate", cxTwoPointCopy)
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
    toolbox.register("select", tools.selTournament, tournsize=3)
    
    random.seed(64)
    
    pop = toolbox.population(n=int(Population_Nbr))
    
    # Numpy equality function (operators.eq) between two arrays returns the
    # equality element wise, which raises an exception in the if similar()
    # check of the hall of fame. Using a different equality function like
    # numpy.array_equal or numpy.allclose solve this issue.
    hof = tools.HallOfFame(1, similar=numpy.array_equal)
    
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", numpy.mean)
    stats.register("std", numpy.std)
    stats.register("min", numpy.min)
    stats.register("max", numpy.max)
    
    algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=int(Max_gen), stats=stats,
                        halloffame=hof)

    return pop, stats, hof

if __name__ == "__main__":
    Nbr_ploy = 150    
    Poly_position_art = np.random.choice(range(2900),size=Nbr_ploy,replace=False)
    sum_signal_art = sumSignal()
    GA()
    
    
    
    
    
    