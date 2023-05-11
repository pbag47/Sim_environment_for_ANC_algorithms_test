close all
clear variables

Sim_setup.Real_filter_sample_rate = 22050 ; % Hz
Sim_setup.ANC_sample_rate = 3000 ; % Hz

%% Searching for minimum convergence time
% Parameters.RLS.lambda = [1.0004] ;
% Parameters.RLS_v2.lambda = [0.737] ;
Parameters.TDLMS.lambda = [0.9924] ;

%% Searching for minimum residuals
% Parameters.RLS.lambda = [0.9988] ;
% Parameters.RLS_v2.lambda = [0.994] ;

%% Tests
% Parameters.RLS.lambda = [1, 0.9988] ;   % 1 | 0.9988
% Parameters.RLS_v2.lambda = [0.737, 0.994] ;  % 0.737 | 0.994
% Parameters.TDLMS.lambda = [0.9924, 0.9983] ;  % 0.9924 | 0.9983
% Parameters.TDLMS_v2.lambda = [0.998] ; % 0.985 | 

%% Sweeping values for lambda
% Parameters.RLS.lambda = linspace(0.9984, 1.0004, 20) ; % min 0.99835 | max 1.0004
% Parameters.RLS_v2.lambda = linspace(0.737, 0.994, 20) ; % min 0.737 | max 0.994
% Parameters.TDLMS.lambda = linspace(0.983, 1-(1e-16), 20) ; % min 0.983 | max 1-(1e-16)
% Parameters.TDLMS_v2.lambda = linspace(0.603, 0.998, 20) ; % min 0.603 | max 0.998 
% Parameters.TDLMS_v3.lambda = linspace(0.885, 0.97, 20) ; % min 0.888 | max 0.9696
% Parameters.LMS.mu = linspace(0.5, 6, 50) ;

%% Save mode
% Choose to either store simulation results in an external file, or to
% discard the results at the end of the program
save_mode = false ; % Boolean, default: true
if save_mode
    data_file = 'simulation_results_v2.mat' ;
    Parameters = Parse_existing_results(data_file,...
        Sim_setup, Parameters) ;
else
    plot_all_error_curves = false ;
end

%% Input signal and real filter FIR import + conditioning
% Import the UAV audio data @ 48 kHz
load('audio_data.mat')
% Import "real" FIR filter @ 22,050 kHz
load('real_filter.mat')

% Resample UAV audio data from 48 kHz to 22,050 kHz
UAV_audio_data = resample(UAV_noise.Mini.mono,...
    Sim_setup.Real_filter_sample_rate, fs) ;

% UAV_audio_data = randn(size(UAV_audio_data)) ;
% UAV_audio_data = sin(2*pi*440 * 1:length(UAV_audio_data) / Sim_setup.Real_filter_sample_rate)' ;
% UAV_audio_data = randn(size(UAV_audio_data)) .* sin(2*pi*440 / ...
%     Sim_setup.Real_filter_sample_rate * 1:length(UAV_audio_data))' ;


%% Desired signal acquisition @ 22,050 kHz ("Real filter" sample rate)
% Filter the resampled audio by the 'real filter'
Filtered_UAV_audio = zeros(length(UAV_audio_data), 1) ;
Buffer = zeros(1, length(Sh)) ;
for i = 1:length(UAV_audio_data)
    Buffer = [UAV_audio_data(i) Buffer(1:length(Buffer)-1)] ;
    Filtered_UAV_audio(i) = Buffer * Sh ;
end

%% Signal conditionning
% Downsample UAV audio data from 22.050 kHz to 3 kHz
Input = resample(UAV_audio_data, Sim_setup.ANC_sample_rate,...
    Sim_setup.Real_filter_sample_rate) ;

% Downsample Filtered_UAV_audio from 22.050 kHz to 3 kHz
Expected_output = resample(Filtered_UAV_audio, Sim_setup.ANC_sample_rate,...
    Sim_setup.Real_filter_sample_rate) ;

%% Assign the start sample from which the algorithm is allowed to work
% The start sample is chosen empirically according to the input signal
% power variations: the aim is to start the algorithms when the power of
% the input signal is relatively constant. 
% If the input power decreases, it can be misinterpreted as an algorihm convergence
% If the input power increases, it can be misinterpreted as an algorithm divergence
ANC_start_sample = 501 + round(length(Sh) * Sim_setup.ANC_sample_rate / Sim_setup.Real_filter_sample_rate) ;

%% Audio and visual renderings of the input signal and the expected output
% figure(500)
% subplot(2,1,1)
% plot(UAV_audio_data)
% subplot(2,1,2)
% plot(Filtered_UAV_audio)

% disp('Sound')
% sound(UAV_noise.Mavic_air_2.mono, fs)
% pause(10)

% sound(UAV_audio_data, 22050)
% pause(10)
% sound(Filtered_UAV_audio, 22050)

% figure(501)
% hold on
% plot(Input)
% hold on
% plot(Expected_output)

%% Algorithm tests
if save_mode
    Results = Algorithm_test(Input, Expected_output, ANC_start_sample, Parameters) ;
    Save_results(data_file, Sim_setup, Results)
    Print_results(data_file)
else
    Print_specific_case_results(Input, Expected_output, ANC_start_sample, Parameters, Sim_setup, plot_all_error_curves) ;
end

