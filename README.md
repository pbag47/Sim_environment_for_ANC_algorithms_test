# Sim_environment_for_ANC_algorithms_test
Matlab code to test Active Noise Cancellation algorithms.

## Objectives

The aim of this simulation environment is to evaluate the performance of several Active Noise Cancellation algorithms to attenuate the noise emitted by a UAV.
Efforts have been made to allow an easy extension of this environment, for example to quickly integrate new ANC algorithms or to change the reference UAV noise.

## Requirements

This program has been developped using Matlab R2017b. It should also be compatible with newer versions of Matlab.

## Program description



## Adding your own ANC algorithm

The next procedure indicates how to add a new algorithm to the simulation environment:

### Algorithm function file
Create a new Matlab function file that contains the new algorithm. 
The name of this file has to be set using the following rule:

Let [Name] denote for your new algorithm name that has been chosen with respects to the file-naming and function-naming formats of Matlab (no space, no special characters except _, valid number of characters).
The function name shall then be:

[Name]_algorithm

Then, save this function file as [Name]_algorithm.m in the same folder than this repository.

### Function inputs
The function inputs are the same for every algorithms, they must follow the next ordered list:
- Input : input signal (column vector)
- Expected_result : output of the "real" FIR filter (column vector)
- ANC_start_sample : sample number at which the ANC system is switched on (int)
- filter_length : adaptive FIR filter length (int)
- variables : algorithm-specific variables (column vector)

### Function output
The function output is the same for every algorithms, it is a "meta-variable" that contains:
- Error : the error signal (column vector)
- t : the computing time (scalar)

### Example
To summarize the previous rules, let's take an example: the implementation of the Least Mean Squares (LMS) algorithm.

In this case, [Name] is set to LMS, which leads to the following function name: LMS_algorithm. 
This function is saved in a Matlab file called LMS_algorithm.m which is displayed as follows.

```Matlab
function [Error, t] = LMS_algorithm(Input, Expected_result, ANC_start_sample, filter_length, variables)
  % :
  % LMS algorithm code
  % :
end
```






