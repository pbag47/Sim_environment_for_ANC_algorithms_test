# Sim_environment_for_ANC_algorithms_test
Matlab code to test Active Noise Cancellation algorithms.

## Objectives

The aim of this simulation environment is to evaluate the performance of several Active Noise Cancellation algorithms to attenuate the noise emitted by a UAV.
Efforts have been made to allow an easy extension of this environment, for example to quickly integrate new ANC algorithms or to change the reference UAV noise.

## Requirements

This program has been developped using standard Matlab R2017b (no additionnal toolbox required). It should also be compatible with newer versions of Matlab.

## Program description



## Adding your own ANC algorithm

The next procedure indicates how to add a new algorithm to the simulation environment:

### Algorithm function file
Create a new Matlab function file that contains the new algorithm. 
The name of this file has to be set using the following rule:

Let ```[Name]``` denote for your new algorithm name that has been chosen with respects to the file-naming and function-naming formats of Matlab (no space, no special characters except _, valid number of characters).
The function name shall then be:

```[Name]```_algorithm

Then, save this function file as ```[Name]```_algorithm.m in the same folder than this repository.

### Function inputs
The function inputs are the same for every algorithms, they must follow the next ordered list:
- ```Input``` : input signal (column vector)
- ```Expected_result``` : output of the "real" FIR filter (column vector)
- ```ANC_start_sample``` : sample number at which the ANC system is switched on (int)
- ```filter_length``` : adaptive FIR filter length (int)
- ```variables``` : algorithm-specific variables (column vector)

### Function output
The function output is the same for every algorithms, it is a "meta-variable" that contains:
- ```Error``` : the error signal (column vector)
- ```t``` : the computing time, in seconds (scalar)

### Example
To summarize the previous rules, let's take an example: the implementation of the Least Mean Squares (LMS) algorithm.

In this case, [Name] is set to LMS, which leads to the following function name: LMS_algorithm. 
This function is saved in a Matlab file called LMS_algorithm.m which is displayed as follows.

```Matlab
function [Error, t] = LMS_algorithm(Input, Expected_result, ANC_start_sample, filter_length, variables)
  
    for i = 1:length(Input)
        % :
        % LMS algorithm code
        % :
    end
  
end
```

The LMS algorithm requires the use of several variables:
- ```mu``` is the convergence factor (scalar), it is an algorithm-specific variable that can be retrieved from the simulation settings using the ```variables``` function input
- ```X``` is a buffer containing the most recent samples of input signal (line-vector)
- ```H``` is the impulse response of the adaptive FIR filter that is controlled by the LMS algorithm (column-vector)

Its implementation leads to:

```Matlab
function [Error, t] = LMS_algorithm(Input, Expected_result, ANC_start_sample, filter_length, variables)
    mu = variables(1) ;
    Error = zeros(length(Input), 1) ;
    X = zeros(1, filter_length) ;
    H = zeros(filter_length, 1) ;
  
    for i = 1:length(Input)
        X = [Input(i) X(1 : filter_length-1] ;
        Error(i) = Expected_result(i) - X*H ;
        H = H + mu * Error(i) * X' ;
    end
end
```

### Additionnal features
#### Computation time
The simulation environment allows to compare the computing time of the tested algorithms. 
The function output ```t``` is dedicated to store the computing time of the algorithm, using ```tic()``` and ```toc()``` functions.

Using the previously-defined example file, this feature can be implemented as follows:

```Matlab
function [Error, t] = LMS_algorithm(Input, Expected_result, ANC_start_sample, filter_length, variables)
    mu = variables(1) ;
    t = NaN :
    Error = zeros(length(Input), 1) ;
    X = zeros(1, filter_length) ;
    H = zeros(filter_length, 1) ;
    tic()
  
    for i = 1:length(Input)
        X = [Input(i) X(1 : filter_length-1] ;
        Error(i) = Expected_result(i) - X*H ;
        H = H + mu * Error(i) * X' ;
    end
  
    t = toc() ;
    disp(['    Algorithm running time : ', num2str(t), ' s'])
end
```

#### ANC start delay
The algorithm performance evaluation requires to compare the amount of error before and after the algorithm convergence.
To do so, it is necessary to turn off the ANC during the first samples in order to allow the RMS value of the error to be a consistent comparison reference.

```Matlab
function [Error, t] = LMS_algorithm(Input, Expected_result, ANC_start_sample, filter_length, variables)
    mu = variables(1) ;
    t = NaN ;
    Error = zeros(length(Input), 1) ;
    X = zeros(1, filter_length) ;
    H = zeros(filter_length, 1) ;
    tic()
  
    for i = 1:length(Input)
        X = [Input(i) X(1 : filter_length-1] ;
        Error(i) = Expected_result(i) - X*H ;
    
        if i >= ANC_start_sample
            H = H + mu * Error(i) * X' ;
        end
    end
  
    t = toc() ;
    disp(['    Algorithm running time : ', num2str(t), ' s'])
end
```

#### Divergence detection
Since the divergence of the algorithm is considered as a fault, the results in the case of a divergence are not workable and there is no point evaluating every sample.
To speed up the simulation (especially when testing a lot of algorithms and variables combinations), a divergence detector is integrated in the algorithm. It allows to immediately discard the running simulation when the algorithm diverges, and proceed to the next simulation as soon as possible.

This divergence detector is implemented as follows:

```Matlab
function [Error, t] = LMS_algorithm(Input, Expected_result, ANC_start_sample, filter_length, variables)
    mu = variables(1) ;
    t = NaN ;
    Error = zeros(length(Input), 1) ;
    X = zeros(1, filter_length) ;
    H = zeros(filter_length, 1) ;
    tic()
  
    for i = 1:length(Input)
        X = [Input(i) X(1 : filter_length-1] ;
        Error(i) = Expected_result(i) - X*H ;
    
        if i >= ANC_start_sample
            H = H + mu * Error(i) * X' ;
        end
        if isnan(Error(i))
            disp('    Divergence detected')
            return
        end
    end
  
    t = toc() ;
    disp(['    Algorithm running time : ', num2str(t), ' s'])
end
```

### Simulation settings
Once the function file is ready, the declaration of the new algorithm is done in the "main.m" file, using the struct variable ```Parameters```.
First, add a field called ```[Name]``` to this struct variable, and inside this new field ```[Name]```, then add a new field for every single algorithm-specific variable that your algorithm requires.
In these variable-name fields, store all the tested values of variables in a line-vector.

For example, the following line in "main.m" will test the previously-defined LMS_algorithm function with 3 times in a row, using 3 different values for the convergence factor ```mu```. 
```Matlab
Parameters.LMS.mu = [0.01, 0.1, 0.2] ;
```

