function Updated_parameters = Parse_existing_results(data_file_name, Sim_setup, Parameters)
    % Create a struct that has the same fields as Parameters, but does not
    % contain data
    Updated_parameters = Parameters ;
    Algorithms = fieldnames(Parameters) ;
    for i = 1:length(Algorithms)
        Parameter_names = fieldnames(Parameters.(Algorithms{i})) ;
        for j = 1:length(Parameter_names)
            Updated_parameters.(Algorithms{i}).(Parameter_names{j}) = [] ;
        end
    end
    
    % Parse previous simulations results, and fill 'Updated_parameters' 
    % struct only with parameter combinations that do not already exist
    if isfile(data_file_name)
        load(data_file_name, 'Results')
        if isequal(Results.Parameters, Sim_setup)
            for i = 1:length(Algorithms)
                disp(Algorithms{i})
                Parameter_names = fieldnames(Parameters.(Algorithms{i})) ;
                for j = 1:length(Parameters.(Algorithms{i}).(Parameter_names{1}))
                    try
                        Add_combination = 1 ;
                        Candidates_index = 1:length(Results.(Algorithms{i}).(Parameter_names{1})) ;
                        Equal_parameters = zeros(length(Candidates_index), length(Parameter_names)) ;
                        for k = 1:length(Candidates_index)
                            for m = 1:length(Parameter_names)
                                if Results.(Algorithms{i}).(Parameter_names{m})(Candidates_index(k)) == Parameters.(Algorithms{i}).(Parameter_names{m})(j)
                                    Equal_parameters(k, m) = 1 ;
                                end
                            end
                            if all(Equal_parameters(k, :) == 1)
                                disp(['Parameters combination N°', num2str(j),...
                                    ' already found in existing results -> combination removed from new simulations parameters list'])
                                Add_combination = 0 ;
                                break
                            end
                        end
                        if Add_combination == 1
                            disp(['Adding N°', num2str(j), ' ', Algorithms{i}, ' simulation parameters'])
                            for n = 1:length(Parameter_names)
                                Updated_parameters.(Algorithms{i}).(Parameter_names{n}) = [Updated_parameters.(Algorithms{i}).(Parameter_names{n}),...
                                    Parameters.(Algorithms{i}).(Parameter_names{n})(j)] ;
                            end
                        end
                    catch Error
                        if strcmp(Error.identifier, 'MATLAB:nonExistentField')
                            disp(['New algorithm description detected, adding ',...
                                  Algorithms{i}])
                            Updated_parameters.(Algorithms{i}) = Parameters.(Algorithms{i}) ;
                        else
                            throw(Error)
                        end
                    end
                end
            end
        else
            disp(' ---- Warning ---- The current simulation setup differs from the existing simulation results')
            Updated_parameters = Parameters ;
        end
    else
        disp(' ---- Warning ---- No data file found at provided path')
        for i = 1:length(Algorithms)
            Parameter_names = fieldnames(Parameters.(Algorithms{i})) ;
            for j = 1:length(Parameter_names)
                for k = 1:length(Parameters.(Algorithms{i}).(Parameter_names{j}))
                    Updated_parameters.(Algorithms{i}).(Parameter_names{j}) = [Updated_parameters.(Algorithms{i}).(Parameter_names{j}),...
                        Parameters.(Algorithms{i}).(Parameter_names{j})(k)] ;
                end
            end
        end  
        Updated_parameters = Parameters ;
    end
end