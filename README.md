# Scrip_Yang
scrip yang

# files description:
## getSignal.m
* This function is transforming a Start Polymerase position in signal intensities versus time passed.
  * input: 1 polymerase start position
  * output: 1 signal intensity 

## sumSignal.m
* input: a set of polymerase start positions
* output: sum of signals with 400 'frame' whose size is equal to the experimental data. 

## simulation_amp_polyNbr.m
* It is the main program containing:
  * define parameters
  * estimate noise by applying genetic algorithm to 8 experiment data
  * decovolution a artificial data with different amplitude 
  * decovolution a artificial data with different polymerase number
  * local optimization 
  * calculate jaccard distance 

## plot.m 
* Use result from simulation_amp_polyNbr.m to visualize the results.

## ga_iters_test.m
* test how many iterations is the best.

## sumSignal_diff_speed.m
* sumSignal for different polymerase speed. 
## GA_diff_polySpeed_test.m
* use a artificial data to test the performance of different polymerase speed.

## /data_workSpace/
* save results of different simulation.

# What's next
* Genetic algorithm should be written into a function. 
* write the whole procedure for experimental data
