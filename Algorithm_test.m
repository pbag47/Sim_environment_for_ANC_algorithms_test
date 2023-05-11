function Results = Algorithm_test(Input, Expected_output, ANC_start_sample, Parameters)
    Results = struct() ;
    filter_length = 100 ;
    Average_length = 10 ;
    Algorithm_names = fieldnames(Parameters) ;
    curve_number = 0 ;
    header = {} ;
    for i=1:length(Algorithm_names)
        alg_name = Algorithm_names{i} ;
        Variable_names = fieldnames(Parameters.(alg_name)) ;
        for j=1:length(Variable_names)
            Results.(alg_name).(Variable_names{j}) = zeros(1, length(Parameters.(alg_name).(Variable_names{j}))) ;
        end
        
        Results.(alg_name).convergence = zeros(1, length(Parameters.(alg_name).(Variable_names{1}))) ;
        Results.(alg_name).residuals = zeros(1, length(Parameters.(alg_name).(Variable_names{1}))) ;
        Results.(alg_name).computing_time = zeros(1, length(Parameters.(alg_name).(Variable_names{1}))) ;
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
                if residuals > Error_RMS(1) / 4 % 5e-3
                    disp('    Divergence detected')
                    residuals = NaN ;
                    computing_time = NaN ;
                else
                    for k = 1:length(Error_RMS)
                        if Error_RMS(k) < residuals && isnan(convergence)
                            convergence = k+1 ;
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
            
            disp(['    convergence : ', num2str(convergence), 'pts'])
            disp(['    residuals : ', num2str(residuals)])
            
%             figure(1001)
%             hold on
%             plot(Error_RMS)
%             drawnow
%             
%             figure(1000)
%             hold on
%             plot(Error)
%             drawnow
        end
%         figure(1000)
%         hold on
%         legend(header, 'Interpreter', 'none')
%         title('Error vs sample number')
%         xlabel('Sample number')
%         ylabel('Error (V)')
    end
end