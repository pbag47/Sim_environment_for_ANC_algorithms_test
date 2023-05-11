function Print_specific_case_results(Input, Expected_output, ANC_start_sample, Parameters, Sim_setup, plot_all_error_curves)
    filter_length = 100 ;  % 100
    Average_length = 10 ; % Length of the sliding-window RMS value (empirically-chosen)
    
    %% Variables initialization
    Algorithm_names = fieldnames(Parameters) ;
    Results = struct() ;
    curve_number = 0 ;
    header = {} ;
    
    %% Algorithm tests
    for i=1:length(Algorithm_names)
        % Extracts the name of the variables of the currently-tested
        % algorithm and initializes the Results structure
        alg_name = Algorithm_names{i} ;
        Variable_names = fieldnames(Parameters.(alg_name)) ;
        for j=1:length(Variable_names)
            Results.(alg_name).(Variable_names{j}) = zeros(1, length(Parameters.(alg_name).(Variable_names{j}))) ;
        end
        Results.(alg_name).convergence = zeros(1, length(Parameters.(alg_name).(Variable_names{1}))) ;
        Results.(alg_name).residuals = zeros(1, length(Parameters.(alg_name).(Variable_names{1}))) ;
        Results.(alg_name).computing_time = zeros(1, length(Parameters.(alg_name).(Variable_names{1}))) ;
        
        % Runs simulation for each case (each value of the variables)
        for j=1:length(Parameters.(alg_name).(Variable_names{1}))
            curve_number = curve_number + 1 ;
            header{curve_number} = alg_name ;
            var_values = zeros(1, length(Variable_names)) ;
            for k=1:length(Variable_names)
                var_name = Variable_names{k} ;
                header{curve_number} = [header{curve_number}, ', ', var_name, ' = ', num2str(Parameters.(alg_name).(var_name)(j))] ;
                var_values(k) = Parameters.(alg_name).(var_name)(j) ;
            end
            disp(header{curve_number})
            
            function_name = strcat(alg_name, '_algorithm') ;
            algorithm_function = str2func(function_name) ;
            [Error, computing_time] = algorithm_function(Input,...
                Expected_output, ANC_start_sample, filter_length, var_values) ;
            
            Error_RMS = zeros(length(Error)-round(Average_length/2) - ANC_start_sample, 1) ;
            for k = 1:length(Error)-round(Average_length/2) - ANC_start_sample
                if k <= Average_length
                    Error_RMS(k) = mean(rms(Error(ANC_start_sample:ANC_start_sample+Average_length))) ;
                else
                    Error_RMS(k) = mean(rms(Error(k-round(Average_length/2)+ANC_start_sample:k+round(Average_length/2)+ANC_start_sample))) ;
                    if isnan(Error_RMS(k))
                        Error_RMS(k) = 10e50 ;
                    end
                end
            end
            convergence = NaN ;
            residuals = NaN ;
            if sum(Error_RMS > 1) == 0
                number_of_active_ANC_samples = length(Input) - ANC_start_sample ;
                residuals = mean(Error_RMS(end-round(number_of_active_ANC_samples/5):end)) ;
                if residuals > Error_RMS(1) / 10 || residuals == 0 % Error_RMS(1) / 4 % 5e-3
                    disp('    Divergence or too slow convrgence detected')
                    residuals = NaN ;
                    computing_time = NaN ;
                else
                    for k = 1:length(Error_RMS)
                        if Error_RMS(k) < residuals && isnan(convergence)
                            convergence = k ;
                        end
%                     for k = length(Error_RMS):-1:1
%                         if Error_RMS(k) > (Error_RMS(1)+residuals)/4 && isnan(convergence)
%                             convergence = k+1 ;
%                         end
                        if Error_RMS(k) > 1.5 * Error_RMS(1)
                            disp('    Divergence detected')
                            convergence = NaN ;
                            residuals = NaN ;
                            computing_time = NaN ;
                            break
                        end
                    end
                end
            else
                disp('    Divergence detected: error RMS reached a too high value')
                computing_time = NaN ;
            end
            for k=1:length(Variable_names)
                Results.(alg_name).(Variable_names{k})(j) = Parameters.(alg_name).(Variable_names{k})(j) ;
            end
            Results.(alg_name).convergence(j) = convergence ;
            Results.(alg_name).residuals(j) = residuals ;
            Results.(alg_name).computing_time(j) = computing_time ;
            
            %% Display run-specific results
            
            disp(['    convergence : ', num2str(convergence), 'pts'])
            disp(['    residuals : ', num2str(residuals)])
            
            if plot_all_error_curves
                figure(1000)
                hold on
                plot(Error)

                figure(1001)
                hold on
                plot(Error_RMS)
            end
        end
    end
    curve_number = curve_number + 1 ;
    header{curve_number} = 'ANC start sample' ;
    
    if plot_all_error_curves
        figure(1000)
        hold on
        plot([ANC_start_sample, ANC_start_sample], [min(Expected_output), max(Expected_output)], '--k')
        legend(header, 'Interpreter', 'none')
        title('Error vs sample number')
        xlabel('Sample number')
        ylabel('Error (V)')

        figure(1001)
        hold on
        plot([Average_length, Average_length], [0, rms(Expected_output)], '--k')
        legend(header, 'Interpreter', 'none')
        title('Error RMS vs RMS sample number')
        xlabel('RMS sample number')
        ylabel('Error RMS (V_R_M_S)')
    end
    
    
    %% Display results summary
    
    % Extracts the names of tested algorithms from the results structure
    Results_fieldnames = fieldnames(Results) ;
    Algorithm_names = {} ;
    for i=1:length(Results_fieldnames)
        if ~strcmp(Results_fieldnames{i}, 'Parameters')
            Algorithm_names{length(Algorithm_names)+1, 1} = Results_fieldnames{i} ;
        end
    end
    
    % Formats the simulation results into 5 variables for boxplots:
    %   -   "full_cv" concatenates the convergence results series of every 
    %       algorithm into a single series
    %   -   "full_res" concatenates the residuals results series of every 
    %       algorithm into a single series
    %   -   "full_ct" concatenates the computing time results series of every 
    %       algorithm into a single series
    %   -   "full_names" concatenates the algorithm name series of every 
    %       algorithm into a single series
    %   -   "full_x_label_coordinates" is used to scatter-plot every 
    %       simulation result as a cross on the boxplot. To do so, an 
    %       equivalence is made between each algorithm name and its 
    %       respective x-axis coordinate on the boxplot graph.
    convergence = zeros(1, length(Algorithm_names)) ;
    residuals = zeros(1, length(Algorithm_names)) ;
    computing_time = zeros(1, length(Algorithm_names)) ;
    full_cv = [] ;
    full_res = [] ;
    full_ct = [] ;
    full_names = [] ;
    full_x_label_coordinates = [] ;
    for i=1:length(Algorithm_names)
        cv = (Results.(Algorithm_names{i}).convergence/Sim_setup.ANC_sample_rate)' ;
        res = (Results.(Algorithm_names{i}).residuals)' ;
        ct = (Results.(Algorithm_names{i}).computing_time)' ;
        names = repmat({Algorithm_names{i}}, [length(cv), 1]) ;
        x_label_coordinates = repmat(i, [length(cv), 1]) ;
        convergence(i) = min(Results.(Algorithm_names{i}).convergence)/Sim_setup.ANC_sample_rate ;
        residuals(i) = min(Results.(Algorithm_names{i}).residuals) ;
        computing_time(i) = min(Results.(Algorithm_names{i}).computing_time) ;
        full_cv = [full_cv ; cv] ;
        full_res = [full_res ; res] ;
        full_ct = [full_ct ; ct] ;
        full_names = [full_names ; names] ;
        full_x_label_coordinates = [full_x_label_coordinates ; x_label_coordinates] ;
    end
    
    figure(1)
    subplot(3,1,1)
    boxplot(full_cv, full_names, 'Whisker', Inf)
    hold on
    scatter(full_x_label_coordinates, full_cv, 20, 'x', 'jitter', 'on', 'jitterAmount', 0.05)
    ylim([0, Inf])
    title('Convergence time comparison')
    ylabel('Convergence time (s)')
    set(gca, 'TickLabelInterpreter', 'none')
    
    figure(1)
    subplot(3,1,2)
    boxplot(full_res, full_names, 'Whisker', Inf)
    hold on
    scatter(full_x_label_coordinates, full_res, 20, 'x', 'jitter', 'on', 'jitterAmount', 0.05)
    ylim([0, Inf])
    title('Residuals comparison')
    ylabel('Residuals (V_R_M_S)')
    set(gca, 'TickLabelInterpreter', 'none')
    
    figure(1)
    subplot(3,1,3)
    boxplot(full_ct, full_names, 'Whisker', Inf)
    hold on
    scatter(full_x_label_coordinates, full_ct, 20, 'x', 'jitter', 'on', 'jitterAmount', 0.05)
    ylim([0, Inf])
    title('Computing time comparison')
    ylabel('Computing time (s)')
    set(gca, 'TickLabelInterpreter', 'none')
    
    % For each algorithm, plots the convergence, residuals and computing
    % time results vs a tested parameter (depends on the algorithm)
    for i=1:length(Algorithm_names)
        % Extracts the tested parameters names from the results structure
        Variable_names = {} ;
        algorithm_fieldnames = fieldnames(Results.(Algorithm_names{i})) ;
        for j=1:length(algorithm_fieldnames)
            if ~strcmp(algorithm_fieldnames{j}, 'convergence') && ...
                    ~strcmp(algorithm_fieldnames{j}, 'residuals') && ...
                    ~strcmp(algorithm_fieldnames{j}, 'computing_time')
                Variable_names{length(Variable_names)+1, 1} = algorithm_fieldnames{j} ;
            end
        end
        
        figure(i+1)
        % This program only works for a single-variable sweep since the 
        % results visualization of a n-variables sweep case requires 
        % graphs with at least n+1 dimensions.
        if length(Variable_names) == 1
            % Finds the minimum value of each series to highlight its
            % position with a marker
            [~, cv_min_index] = min(Results.(Algorithm_names{i}).convergence/Sim_setup.ANC_sample_rate) ;
            [~, r_min_index] = min(Results.(Algorithm_names{i}).residuals) ;
            [~, ct_min_index] = min(Results.(Algorithm_names{i}).computing_time) ;
            
            % Convergence time series
            subplot(3,1,1)
            plot(Results.(Algorithm_names{i}).(Variable_names{1}), ...
                (Results.(Algorithm_names{i}).convergence/Sim_setup.ANC_sample_rate),...
                'LineStyle', '-', 'Marker', '.')
            hold on
            scatter(Results.(Algorithm_names{i}).(Variable_names{1})(cv_min_index), ...
                Results.(Algorithm_names{i}).convergence(cv_min_index)/Sim_setup.ANC_sample_rate)
            ylim([0, Inf])
            title([Algorithm_names{i}, ' convergence time vs ', Variable_names{1}], 'Interpreter', 'none')
            xlabel(Variable_names{1}, 'Interpreter', 'none')
            ylabel('Convergence time (s)')
            
            % Residuals series
            subplot(3,1,2)
            plot(Results.(Algorithm_names{i}).(Variable_names{1}), ...
                Results.(Algorithm_names{i}).residuals,...
                'LineStyle', '-', 'Marker', '.')
            hold on
            scatter(Results.(Algorithm_names{i}).(Variable_names{1})(r_min_index), ...
                Results.(Algorithm_names{i}).residuals(r_min_index))
            ylim([0, Inf])
            title([Algorithm_names{i}, ' residuals vs ', Variable_names{1}], 'Interpreter', 'none')
            xlabel(Variable_names{1}, 'Interpreter', 'none')
            ylabel('Residuals (V_R_M_S)')
            
            % Computing time series
            subplot(3,1,3)
            plot(Results.(Algorithm_names{i}).(Variable_names{1}), ...
                Results.(Algorithm_names{i}).computing_time,...
                'LineStyle', '-', 'Marker', '.')
            hold on
            scatter(Results.(Algorithm_names{i}).(Variable_names{1})(ct_min_index), ...
                Results.(Algorithm_names{i}).computing_time(ct_min_index))
            ylim([0, Inf])
            title([Algorithm_names{i}, ' computing time vs ', Variable_names{1}], 'Interpreter', 'none')
            xlabel(Variable_names{1}, 'Interpreter', 'none')
            ylabel('Computing time (s)')
        else
            disp(['---- Warning ----, ', Algorithm_names{i}, ...
                ' invalid number of variables to display detailed results graph'])
        end
    end
    
    if plot_all_error_curves
        figure(1000)
        plotbrowser
        figure(1001)
        plotbrowser
    end
end