function Save_results(data_file_name, Sim_setup, Results_to_save)
%     global Resultats
    if isfile(data_file_name)
        load(data_file_name, 'Results')
        if isequal(Results.Parameters, Sim_setup)
            Evaluated_algorithms = fieldnames(Results_to_save) ;
            for i = 1:length(Evaluated_algorithms)
                if ~strcmp(Evaluated_algorithms{i}, 'Parameters')
                    try
                        Var_names = fieldnames(Results.(Evaluated_algorithms{i})) ;
                    catch
                        % Allows the creation of storage locations for
                        % newly encountered variables, if simulation
                        % results are available
                        Var_names = fieldnames(Results_to_save.(Evaluated_algorithms{i})) ;
                        if any(strcmp(Var_names, 'convergence'))
                            for j = 1:length(Var_names)
                                for k = 1:length(Results_to_save.(Evaluated_algorithms{i}).(Var_names{j}))
                                    Results.(Evaluated_algorithms{i}).(Var_names{j}) = [] ;
                                end
                            end
                        end
                    end
                    for j = 1:length(Var_names)
                        for k = 1:length(Results_to_save.(Evaluated_algorithms{i}).(Var_names{j}))
                            Results.(Evaluated_algorithms{i}).(Var_names{j}) = [Results.(Evaluated_algorithms{i}).(Var_names{j}),...
                                Results_to_save.(Evaluated_algorithms{i}).(Var_names{j})(k)] ;
                        end
                    end
                end
            end
            % Sort the results to a parameter-ascending order
            Evaluated_algorithms = fieldnames(Results) ;
            for i = 1:length(Evaluated_algorithms)
                if ~strcmp(Evaluated_algorithms(i), 'Parameters')
                    Parameter_names = fieldnames(Results.(Evaluated_algorithms{i})) ;
                    [Results.(Evaluated_algorithms{i}).(Parameter_names{1}), I] = sort(Results.(Evaluated_algorithms{i}).(Parameter_names{1})) ;
                    for j = 2:length(Parameter_names)
                        Results.(Evaluated_algorithms{i}).(Parameter_names{j}) = Results.(Evaluated_algorithms{i}).(Parameter_names{j})(I) ;
                    end
                end
            end
%             Resultats = Results ;
            save(data_file_name, 'Results')
            disp('Data saved')
        else
            disp(' ---- Warning ---- No existing results with matching parameters',...
                'found at provided file path, creating a new file to save data...')
            Results_to_save.Parameters = Sim_setup ;
            Results = Results_to_save ;
%             Resultats = Results ;
            save('data_recovery_file.mat', 'Results')
            disp('Results saved to data_recovery_file.mat')
        end
    else
        disp(' ---- Warning ---- No data file found at provided path, creating a new file to save data...')
        Results_to_save.Parameters = Sim_setup ;
        Results = Results_to_save ;
        save(data_file_name, 'Results')
        disp('File created, data saved')
    end
end