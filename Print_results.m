function Print_results(filename)
    load(filename, 'Results')
    
    Results_fieldnames = fieldnames(Results) ;
    Algorithm_names = {} ;
    for i=1:length(Results_fieldnames)
        if ~strcmp(Results_fieldnames{i}, 'Parameters')
            Algorithm_names{length(Algorithm_names)+1, 1} = Results_fieldnames{i} ;
        end
    end
    
    X = categorical(Algorithm_names) ;
    X = reordercats(X, Algorithm_names) ;
    convergence = zeros(1, length(Algorithm_names)) ;
    residuals = zeros(1, length(Algorithm_names)) ;
    computing_time = zeros(1, length(Algorithm_names)) ;
    for i=1:length(Algorithm_names)
        convergence(i) = min(Results.(Algorithm_names{i}).convergence)/Results.Parameters.ANC_sample_rate ;
        residuals(i) = min(Results.(Algorithm_names{i}).residuals) ;
        computing_time(i) = min(Results.(Algorithm_names{i}).computing_time) ;
    end
    
    figure(1)
    subplot(3,1,1)
    bar(X, convergence)
    ylim([0, Inf])
    title('Minimum convergence time comparison')
    ylabel('Convergence time (s)')
    set(gca, 'TickLabelInterpreter', 'none')
    
    subplot(3,1,2)
    bar(X, residuals)
    ylim([0, Inf])
    title('Minimum residuals comparison')
    ylabel('Residuals (V_R_M_S)')
    set(gca, 'TickLabelInterpreter', 'none')
    
    subplot(3,1,3)
    bar(X, computing_time)
    ylim([0, Inf])
    title('Minimum computing time comparison')
    ylabel('Computing time (s)')
    set(gca, 'TickLabelInterpreter', 'none')
    
    for i=1:length(Algorithm_names)
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
        if length(Variable_names) == 1
            [~, cv_min_index] = min(Results.(Algorithm_names{i}).convergence/Results.Parameters.ANC_sample_rate) ;
            [~, r_min_index] = min(Results.(Algorithm_names{i}).residuals) ;
            [~, ct_min_index] = min(Results.(Algorithm_names{i}).computing_time) ;
            subplot(3,1,1)
            plot(Results.(Algorithm_names{i}).(Variable_names{1}), ...
                (Results.(Algorithm_names{i}).convergence/Results.Parameters.ANC_sample_rate),...
                'LineStyle', '-', 'Marker', '.')
            hold on
            scatter(Results.(Algorithm_names{i}).(Variable_names{1})(cv_min_index), ...
                Results.(Algorithm_names{i}).convergence(cv_min_index)/Results.Parameters.ANC_sample_rate)
            ylim([0, Inf])
            title([Algorithm_names{i}, ' convergence time vs ', Variable_names{1}], 'Interpreter', 'none')
            xlabel(Variable_names{1}, 'Interpreter', 'none')
            ylabel('Convergence time (s)')

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
    
    
    
    
    
    %% 1st version
%     
%     figure(1)
%     subplot(2,1,1)
%     X = categorical({'LMS', 'NLMS', 'HLMS', 'Leaky-LMS', 'Leaky-NLMS', 'LMP'});
%     X = reordercats(X,{'LMS', 'NLMS', 'HLMS', 'Leaky-LMS', 'Leaky-NLMS', 'LMP'});
%     Y = [min(Results.LMS.convergence)/Results.Parameters.ANC_sample_rate,...
%         min(Results.NLMS.convergence)/Results.Parameters.ANC_sample_rate,...
%         min(Results.HLMS.convergence)/Results.Parameters.ANC_sample_rate,...
%         min(Results.LLMS.convergence)/Results.Parameters.ANC_sample_rate,...
%         min(Results.LNLMS.convergence)/Results.Parameters.ANC_sample_rate,...
%         min(Results.LMP.convergence)/Results.Parameters.ANC_sample_rate] ;
%     bar(X, Y)
%     ylim([0, Inf])
%     title('Minimum convergence time comparison')
%     ylabel('Convergence time (s)')
%     
%     subplot(2,1,2)
%     X = categorical({'LMS', 'NLMS', 'HLMS', 'Leaky-LMS', 'Leaky-NLMS', 'LMP'});
%     X = reordercats(X,{'LMS', 'NLMS', 'HLMS', 'Leaky-LMS', 'Leaky-NLMS', 'LMP'});
%     Y = [min(Results.LMS.residuals),...
%         min(Results.NLMS.residuals),...
%         min(Results.HLMS.residuals),...
%         min(Results.LLMS.residuals),...
%         min(Results.LNLMS.residuals),...
%         min(Results.LMP.residuals)] ;
%     bar(X, Y)
%     ylim([0, Inf])
%     title('Minimum residuals comparison')
%     ylabel('Residuals (V_R_M_S)')
%     
%     figure(2)
%     [~, c_min_index] = min(Results.LMS.convergence/Results.Parameters.ANC_sample_rate) ;
%     [~, r_min_index] = min(Results.LMS.residuals) ;
%     subplot(2,1,1)
%     plot(Results.LMS.mu, (Results.LMS.convergence/Results.Parameters.ANC_sample_rate),...
%         'LineStyle', '-', 'Marker', '.')
%     hold on
%     scatter(Results.LMS.mu(c_min_index), Results.LMS.convergence(c_min_index)/Results.Parameters.ANC_sample_rate)
%     ylim([0, Inf])
%     title('LMS convergence time vs \mu')
%     xlabel('\mu')
%     ylabel('Convergence time (s)')
%     legend('Convergence time vs \mu', 'Minimum convergence time', 'Minimum residuals')
%     subplot(2,1,2)
%     plot(Results.LMS.mu, Results.LMS.residuals,...
%         'LineStyle', '-', 'Marker', '.')
%     hold on
%     scatter(Results.LMS.mu(r_min_index), Results.LMS.residuals(r_min_index))
%     ylim([0, Inf])
%     title('LMS residuals vs \mu')
%     xlabel('\mu')
%     ylabel('Residuals (V_R_M_S)')
%     
%     figure(3)
%     [~, c_min_index] = min(Results.NLMS.convergence/Results.Parameters.ANC_sample_rate) ;
%     [~, r_min_index] = min(Results.NLMS.residuals) ;
%     subplot(2,1,1)
%     plot(Results.NLMS.mu, (Results.NLMS.convergence/Results.Parameters.ANC_sample_rate),...
%         'LineStyle', '-', 'Marker', '.')
%     hold on
%     scatter(Results.NLMS.mu(c_min_index), Results.NLMS.convergence(c_min_index)/Results.Parameters.ANC_sample_rate)
%     ylim([0, Inf])
%     title('NLMS convergence time vs \mu')
%     xlabel('\mu')
%     ylabel('Convergence time (s)')
%     subplot(2,1,2)
%     plot(Results.NLMS.mu, Results.NLMS.residuals,...
%         'LineStyle', '-', 'Marker', '.')
%     hold on
%     scatter(Results.NLMS.mu(r_min_index), Results.NLMS.residuals(r_min_index))
%     ylim([0, Inf])
%     title('NLMS residuals vs \mu')
%     xlabel('\mu')
%     ylabel('Residuals (V_R_M_S)')
%     
% %     figure(4) 
% %     subplot(2,1,1)
% %     colormap(hsv)
% %     w_series = {} ;
% %     w_series_number = 0 ;
% %     w = [] ;
% %     delta = [] ;
% %     delta_series = {} ;
% %     delta_series_number = 0 ;
% %     for j = 1:length(Results.RLS.w)
% %         if ~any(w == Results.RLS.w(j))
% %             I = find(Results.RLS.w == Results.RLS.w(j)) ;
% %             w = [w, Results.RLS.w(j)] ;
% %             w_series_number = w_series_number + 1 ;
% %             w_s = [Results.RLS.delta(I) ; 
% %                 Results.HLMS.convergence(I)/Results.Parameters.ANC_sample_rate ;
% %                 Results.HLMS.residuals(I)] ;
% %             [w_s(1,:), I2] = sort(w_s(1, :)) ;
% %             w_s(2, :) = w_s(2, I2) ;
% %             w_s(3,:) = w_s(3, I2) ;
% %             w_series{w_series_number, 1} = [nan, w_s(1,:), nan ;
% %                 nan, w_s(2,:), nan ;
% %                 nan, w_s(3,:), nan] ;
% %             hold on
% %             w_axis = Results.RLS.w(j) * ones(1, length(I)+2) ;
% %             w_axis(1) = nan ;
% %             w_axis(length(I)+2) = nan ;
% %             color_values = w_series{w_series_number, 1}(2, :) ;
% %             for k = 1:length(color_values)
% %                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
% %                     color_values(k) = 1 ;
% %                 end
% %             end
% %             patch(w_axis, w_series{w_series_number, 1}(1, :), w_series{w_series_number, 1}(2, :),...
% %                 color_values,...
% %                 'FaceColor', 'none', 'EdgeColor', 'interp')
% %             colorbar
% %             view(3)
% %         end
% %     end
% %     for j = 1:length(Results.RLS.delta)
% %         if ~any(delta == Results.RLS.delta(j))
% %             I = find(Results.RLS.delta == Results.RLS.delta(j)) ;
% %             delta = [delta, Results.RLS.delta(j)] ;
% %             delta_series_number = delta_series_number + 1 ;
% %             delta_s = [Results.RLS.w(I) ; 
% %                 Results.HLMS.convergence(I)/Results.Parameters.ANC_sample_rate ;
% %                 Results.HLMS.residuals(I)] ;
% %             [delta_s(1,:), I2] = sort(delta_s(1, :)) ;
% %             delta_s(2, :) = delta_s(2, I2) ;
% %             delta_s(3,:) = delta_s(3, I2) ;
% %             delta_series{delta_series_number, 1} = [nan, delta_s(1,:), nan ;
% %                 nan, delta_s(2,:), nan ;
% %                 nan, delta_s(3,:), nan] ;
% %             hold on
% %             delta_axis = Results.RLS.delta(j) * ones(1, length(I)+2) ;
% %             delta_axis(1) = nan ;
% %             delta_axis(length(I)+2) = nan ;
% %             color_values = delta_series{delta_series_number, 1}(2, :) ;
% %             for k = 1:length(color_values)
% %                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
% %                     color_values(k) = 1 ;
% %                 end
% %             end
% %             patch(delta_series{delta_series_number, 1}(1, :), delta_axis, delta_series{delta_series_number, 1}(2, :),...
% %                 color_values,...
% %                 'FaceColor', 'none', 'EdgeColor', 'interp')
% %             colorbar
% %             view(3)
% %         end
% %     end
% %     hold on
% %     scatter3(Results.RLS.w, Results.RLS.delta, Results.RLS.convergence/Results.Parameters.ANC_sample_rate)
% %     xlabel('w')
% %     ylabel('\delta')
% %     zlabel('Convergence time (s)')
% %     title('RLS convergence time vs w and \delta')
% %     subplot(2,1,2)
% %     stem3(Results.RLS.w, Results.RLS.delta, Results.RLS.residuals)
% %     xlabel('w')
% %     ylabel('\delta')
% %     zlabel('Residuals (V_R_M_S)')
% %     title('RLS residuals vs w and \delta')
% 
%     figure(5)   
%     mu_series = {} ;
%     mu_series_number = 0 ;
%     mu = [] ;
%     p = [] ;
%     p_series = {} ;
%     p_series_number = 0 ;
%     for j = 1:length(Results.HLMS.mu)
%         if ~any(mu == Results.HLMS.mu(j))
%             I = find(Results.HLMS.mu == Results.HLMS.mu(j)) ;
%             mu = [mu, Results.HLMS.mu(j)] ;
%             mu_series_number = mu_series_number + 1 ;
%             mu_s = [Results.HLMS.p(I) ; 
%                 Results.HLMS.convergence(I)/Results.Parameters.ANC_sample_rate ;
%                 Results.HLMS.residuals(I)] ;
%             [mu_s(1,:), I2] = sort(mu_s(1, :)) ;
%             mu_s(2, :) = mu_s(2, I2) ;
%             mu_s(3,:) = mu_s(3, I2) ;
%             mu_series{mu_series_number, 1} = [nan, mu_s(1,:), nan ;
%                 nan, mu_s(2,:), nan ;
%                 nan, mu_s(3,:), nan] ;
%             mu_axis = Results.HLMS.mu(j) * ones(1, length(I)+2) ;
%             mu_axis(1) = nan ;
%             mu_axis(length(I)+2) = nan ;
%             subplot(2,1,1)
%             hold on
%             color_values = mu_series{mu_series_number, 1}(2, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(mu_axis, mu_series{mu_series_number, 1}(1, :), mu_series{mu_series_number, 1}(2, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%             subplot(2,1,2)
%             hold on
%             color_values = mu_series{mu_series_number, 1}(3, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(mu_axis, mu_series{mu_series_number, 1}(1, :), mu_series{mu_series_number, 1}(3, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%         end
%     end
%     for j = 1:length(Results.HLMS.p)
%         if ~any(p == Results.HLMS.p(j))
%             I = find(Results.HLMS.p == Results.HLMS.p(j)) ;
%             p = [p, Results.HLMS.p(j)] ;
%             p_series_number = p_series_number + 1 ;
%             p_s = [Results.HLMS.mu(I) ; 
%                 Results.HLMS.convergence(I)/Results.Parameters.ANC_sample_rate ;
%                 Results.HLMS.residuals(I)] ;
%             [p_s(1,:), I2] = sort(p_s(1, :)) ;
%             p_s(2, :) = p_s(2, I2) ;
%             p_s(3,:) = p_s(3, I2) ;
%             p_series{p_series_number, 1} = [nan, p_s(1,:), nan ;
%                 nan, p_s(2,:), nan ;
%                 nan, p_s(3,:), nan] ;
%             p_axis = Results.HLMS.p(j) * ones(1, length(I)+2) ;
%             p_axis(1) = nan ;
%             p_axis(length(I)+2) = nan ;
%             subplot(2,1,1)
%             hold on
%             color_values = p_series{p_series_number, 1}(2, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(p_series{p_series_number, 1}(1, :), p_axis, p_series{p_series_number, 1}(2, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%             subplot(2,1,2)
%             hold on
%             color_values = p_series{p_series_number, 1}(3, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(p_series{p_series_number, 1}(1, :), p_axis, p_series{p_series_number, 1}(3, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%         end
%     end
%     subplot(2,1,1)
%     hold on
%     scatter3(Results.HLMS.mu, Results.HLMS.p, Results.HLMS.convergence/Results.Parameters.ANC_sample_rate,...
%         'Marker', '.')
%     colormap(hsv)
%     caxis([min(Results.HLMS.convergence)/Results.Parameters.ANC_sample_rate,...
%         max(Results.HLMS.convergence)/Results.Parameters.ANC_sample_rate])
%     colorbar
%     view(3)
%     xlabel('\mu')
%     ylabel('p')
%     zlabel('Convergence time (s)')
%     title('H-LMS convergence time vs \mu and p')
%     subplot(2,1,2)
%     hold on
%     scatter3(Results.HLMS.mu, Results.HLMS.p, Results.HLMS.residuals,...
%         'Marker', '.')
%     colormap(hsv)
%     caxis([min(Results.HLMS.residuals), max(Results.HLMS.residuals)])
%     colorbar
%     view(3)
%     xlabel('\mu')
%     ylabel('p')
%     zlabel('Residuals (V_R_M_S)')
%     title('H-LMS residuals vs \mu and p')
%     
%     figure(6)   
%     mu_series = {} ;
%     mu_series_number = 0 ;
%     mu = [] ;
%     gamma = [] ;
%     gamma_series = {} ;
%     gamma_series_number = 0 ;
%     for j = 1:length(Results.LLMS.mu)
%         if ~any(mu == Results.LLMS.mu(j))
%             I = find(Results.LLMS.mu == Results.LLMS.mu(j)) ;
%             mu = [mu, Results.LLMS.mu(j)] ;
%             mu_series_number = mu_series_number + 1 ;
%             mu_s = [Results.LLMS.gamma(I) ; 
%                 Results.LLMS.convergence(I)/Results.Parameters.ANC_sample_rate ;
%                 Results.LLMS.residuals(I)] ;
%             [mu_s(1,:), I2] = sort(mu_s(1, :)) ;
%             mu_s(2, :) = mu_s(2, I2) ;
%             mu_s(3,:) = mu_s(3, I2) ;
%             mu_series{mu_series_number, 1} = [nan, mu_s(1,:), nan ;
%                 nan, mu_s(2,:), nan ;
%                 nan, mu_s(3,:), nan] ;
%             mu_axis = Results.LLMS.mu(j) * ones(1, length(I)+2) ;
%             mu_axis(1) = nan ;
%             mu_axis(length(I)+2) = nan ;
%             subplot(2,1,1)
%             hold on
%             color_values = mu_series{mu_series_number, 1}(2, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(mu_axis, mu_series{mu_series_number, 1}(1, :), mu_series{mu_series_number, 1}(2, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%             subplot(2,1,2)
%             hold on
%             color_values = mu_series{mu_series_number, 1}(3, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(mu_axis, mu_series{mu_series_number, 1}(1, :), mu_series{mu_series_number, 1}(3, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%         end
%     end
%     for j = 1:length(Results.LLMS.gamma)
%         if ~any(gamma == Results.LLMS.gamma(j))
%             I = find(Results.LLMS.gamma == Results.LLMS.gamma(j)) ;
%             gamma = [gamma, Results.LLMS.gamma(j)] ;
%             gamma_series_number = gamma_series_number + 1 ;
%             gamma_s = [Results.LLMS.mu(I) ; 
%                 Results.LLMS.convergence(I)/Results.Parameters.ANC_sample_rate ;
%                 Results.LLMS.residuals(I)] ;
%             [gamma_s(1,:), I2] = sort(gamma_s(1, :)) ;
%             gamma_s(2, :) = gamma_s(2, I2) ;
%             gamma_s(3,:) = gamma_s(3, I2) ;
%             gamma_series{gamma_series_number, 1} = [nan, gamma_s(1,:), nan ;
%                 nan, gamma_s(2,:), nan ;
%                 nan, gamma_s(3,:), nan] ;
%             gamma_axis = Results.LLMS.gamma(j) * ones(1, length(I)+2) ;
%             gamma_axis(1) = nan ;
%             gamma_axis(length(I)+2) = nan ;
%             subplot(2,1,1)
%             hold on
%             color_values = gamma_series{gamma_series_number, 1}(2, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(gamma_series{gamma_series_number, 1}(1, :), gamma_axis, gamma_series{gamma_series_number, 1}(2, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%             subplot(2,1,2)
%             hold on
%             color_values = gamma_series{gamma_series_number, 1}(3, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(gamma_series{gamma_series_number, 1}(1, :), gamma_axis, gamma_series{gamma_series_number, 1}(3, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%         end
%     end
%     subplot(2,1,1)
%     hold on
%     scatter3(Results.LLMS.mu, Results.LLMS.gamma, Results.LLMS.convergence/Results.Parameters.ANC_sample_rate,...
%         'Marker', '.')
%     colormap(hsv)
%     caxis([min(Results.HLMS.convergence)/Results.Parameters.ANC_sample_rate,...
%         max(Results.HLMS.convergence)/Results.Parameters.ANC_sample_rate])
%     colorbar
%     view(3)
%     xlabel('\mu')
%     ylabel('\gamma')
%     zlabel('Leaky-LMS convergence time (s)')
%     title('Leaky-LMS convergence time vs \mu and \gamma')
%     subplot(2,1,2)
%     hold on
%     scatter3(Results.LLMS.mu, Results.LLMS.gamma, Results.LLMS.residuals,...
%         'Marker', '.')
%     colormap(hsv)
%     caxis([min(Results.HLMS.residuals), max(Results.HLMS.residuals)])
%     colorbar
%     view(3)
%     xlabel('\mu')
%     ylabel('\gamma')
%     zlabel('Leaky-LMS residuals (V_R_M_S)')
%     title('Leaky-LMS residuals vs \mu and \gamma')
%     
%     figure(7)   
%     mu_series = {} ;
%     mu_series_number = 0 ;
%     mu = [] ;
%     gamma = [] ;
%     gamma_series = {} ;
%     gamma_series_number = 0 ;
%     for j = 1:length(Results.LNLMS.mu)
%         if ~any(mu == Results.LNLMS.mu(j))
%             I = find(Results.LNLMS.mu == Results.LNLMS.mu(j)) ;
%             mu = [mu, Results.LNLMS.mu(j)] ;
%             mu_series_number = mu_series_number + 1 ;
%             mu_s = [Results.LNLMS.gamma(I) ; 
%                 Results.LNLMS.convergence(I)/Results.Parameters.ANC_sample_rate ;
%                 Results.LNLMS.residuals(I)] ;
%             [mu_s(1,:), I2] = sort(mu_s(1, :)) ;
%             mu_s(2, :) = mu_s(2, I2) ;
%             mu_s(3,:) = mu_s(3, I2) ;
%             mu_series{mu_series_number, 1} = [nan, mu_s(1,:), nan ;
%                 nan, mu_s(2,:), nan ;
%                 nan, mu_s(3,:), nan] ;
%             mu_axis = Results.LNLMS.mu(j) * ones(1, length(I)+2) ;
%             mu_axis(1) = nan ;
%             mu_axis(length(I)+2) = nan ;
%             subplot(2,1,1)
%             hold on
%             color_values = mu_series{mu_series_number, 1}(2, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(mu_axis, mu_series{mu_series_number, 1}(1, :), mu_series{mu_series_number, 1}(2, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%             subplot(2,1,2)
%             hold on
%             color_values = mu_series{mu_series_number, 1}(3, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(mu_axis, mu_series{mu_series_number, 1}(1, :), mu_series{mu_series_number, 1}(3, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%         end
%     end
%     for j = 1:length(Results.LNLMS.gamma)
%         if ~any(gamma == Results.LNLMS.gamma(j))
%             I = find(Results.LNLMS.gamma == Results.LNLMS.gamma(j)) ;
%             gamma = [gamma, Results.LNLMS.gamma(j)] ;
%             gamma_series_number = gamma_series_number + 1 ;
%             gamma_s = [Results.LNLMS.mu(I) ; 
%                 Results.LNLMS.convergence(I)/Results.Parameters.ANC_sample_rate ;
%                 Results.LNLMS.residuals(I)] ;
%             [gamma_s(1,:), I2] = sort(gamma_s(1, :)) ;
%             gamma_s(2, :) = gamma_s(2, I2) ;
%             gamma_s(3,:) = gamma_s(3, I2) ;
%             gamma_series{gamma_series_number, 1} = [nan, gamma_s(1,:), nan ;
%                 nan, gamma_s(2,:), nan ;
%                 nan, gamma_s(3,:), nan] ;
%             gamma_axis = Results.LNLMS.gamma(j) * ones(1, length(I)+2) ;
%             gamma_axis(1) = nan ;
%             gamma_axis(length(I)+2) = nan ;
%             subplot(2,1,1)
%             hold on
%             color_values = gamma_series{gamma_series_number, 1}(2, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(gamma_series{gamma_series_number, 1}(1, :), gamma_axis, gamma_series{gamma_series_number, 1}(2, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%             subplot(2,1,2)
%             hold on
%             color_values = gamma_series{gamma_series_number, 1}(3, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(gamma_series{gamma_series_number, 1}(1, :), gamma_axis, gamma_series{gamma_series_number, 1}(3, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%         end
%     end
%     subplot(2,1,1)
%     hold on
%     scatter3(Results.LNLMS.mu, Results.LNLMS.gamma, Results.LNLMS.convergence/Results.Parameters.ANC_sample_rate,...
%         'Marker', '.')
%     colormap(hsv)
%     caxis([min(Results.HLMS.convergence)/Results.Parameters.ANC_sample_rate,...
%         max(Results.HLMS.convergence)/Results.Parameters.ANC_sample_rate])
%     colorbar
%     view(3)
%     xlabel('\mu')
%     ylabel('\gamma')
%     zlabel('Leaky-NLMS convergence time (s)')
%     title('Leaky-NLMS convergence time vs \mu and \gamma')
%     subplot(2,1,2)
%     hold on
%     scatter3(Results.LNLMS.mu, Results.LNLMS.gamma, Results.LNLMS.residuals,...
%         'Marker', '.')
%     colormap(hsv)
%     caxis([min(Results.HLMS.residuals), max(Results.HLMS.residuals)])
%     colorbar
%     view(3)
%     xlabel('\mu')
%     ylabel('\gamma')
%     zlabel('Leaky-NLMS residuals (V_R_M_S)')
%     title('Leaky-NLMS residuals vs \mu and \gamma')
%     
%     figure(8)   
%     mu_series = {} ;
%     mu_series_number = 0 ;
%     mu = [] ;
%     p = [] ;
%     p_series = {} ;
%     p_series_number = 0 ;
%     for j = 1:length(Results.LMP.mu)
%         if ~any(mu == Results.LMP.mu(j))
%             I = find(Results.LMP.mu == Results.LMP.mu(j)) ;
%             mu = [mu, Results.LMP.mu(j)] ;
%             mu_series_number = mu_series_number + 1 ;
%             mu_s = [Results.LMP.p(I) ; 
%                 Results.LMP.convergence(I)/Results.Parameters.ANC_sample_rate ;
%                 Results.LMP.residuals(I)] ;
%             [mu_s(1,:), I2] = sort(mu_s(1, :)) ;
%             mu_s(2, :) = mu_s(2, I2) ;
%             mu_s(3,:) = mu_s(3, I2) ;
%             mu_series{mu_series_number, 1} = [nan, mu_s(1,:), nan ;
%                 nan, mu_s(2,:), nan ;
%                 nan, mu_s(3,:), nan] ;
%             mu_axis = Results.LMP.mu(j) * ones(1, length(I)+2) ;
%             mu_axis(1) = nan ;
%             mu_axis(length(I)+2) = nan ;
%             subplot(2,1,1)
%             hold on
%             color_values = mu_series{mu_series_number, 1}(2, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(mu_axis, mu_series{mu_series_number, 1}(1, :), mu_series{mu_series_number, 1}(2, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%             subplot(2,1,2)
%             hold on
%             color_values = mu_series{mu_series_number, 1}(3, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(mu_axis, mu_series{mu_series_number, 1}(1, :), mu_series{mu_series_number, 1}(3, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%         end
%     end
%     for j = 1:length(Results.LMP.p)
%         if ~any(p == Results.LMP.p(j))
%             I = find(Results.LMP.p == Results.LMP.p(j)) ;
%             p = [p, Results.LMP.p(j)] ;
%             p_series_number = p_series_number + 1 ;
%             p_s = [Results.LMP.mu(I) ; 
%                 Results.LMP.convergence(I)/Results.Parameters.ANC_sample_rate ;
%                 Results.LMP.residuals(I)] ;
%             [p_s(1,:), I2] = sort(p_s(1, :)) ;
%             p_s(2, :) = p_s(2, I2) ;
%             p_s(3,:) = p_s(3, I2) ;
%             p_series{p_series_number, 1} = [nan, p_s(1,:), nan ;
%                 nan, p_s(2,:), nan ;
%                 nan, p_s(3,:), nan] ;
%             p_axis = Results.LMP.p(j) * ones(1, length(I)+2) ;
%             p_axis(1) = nan ;
%             p_axis(length(I)+2) = nan ;
%             subplot(2,1,1)
%             hold on
%             color_values = p_series{p_series_number, 1}(2, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(p_series{p_series_number, 1}(1, :), p_axis, p_series{p_series_number, 1}(2, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%             subplot(2,1,2)
%             hold on
%             color_values = p_series{p_series_number, 1}(3, :) ;
%             for k = 1:length(color_values)
%                 if isnan(color_values(k)) || color_values(k) > 1 || color_values(k) < 0
%                     color_values(k) = 1 ;
%                 end
%             end
%             patch(p_series{p_series_number, 1}(1, :), p_axis, p_series{p_series_number, 1}(3, :),...
%                 color_values,...
%                 'FaceColor', 'none', 'EdgeColor', 'interp')
%         end
%     end
%     subplot(2,1,1)
%     hold on
%     scatter3(Results.LMP.mu, Results.LMP.p, Results.LMP.convergence/Results.Parameters.ANC_sample_rate,...
%         'Marker', '.')
%     colormap(hsv)
%     caxis([min(Results.LMP.convergence)/Results.Parameters.ANC_sample_rate,...
%         max(Results.LMP.convergence)/Results.Parameters.ANC_sample_rate])
%     colorbar
%     view(3)
%     xlabel('\mu')
%     ylabel('p')
%     zlabel('Convergence time (s)')
%     title('LMP convergence time vs \mu and p')
%     subplot(2,1,2)
%     hold on
%     scatter3(Results.LMP.mu, Results.LMP.p, Results.LMP.residuals,...
%         'Marker', '.')
%     colormap(hsv)
%     caxis([min(Results.LMP.residuals), max(Results.LMP.residuals)])
%     colorbar
%     view(3)
%     xlabel('\mu')
%     ylabel('p')
%     zlabel('Residuals (V_R_M_S)')
%     title('LMP residuals vs \mu and p')
end