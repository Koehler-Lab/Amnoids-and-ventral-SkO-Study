classdef TissueStressAnalysis2
    properties (Constant)
        DEFAULT_SMOOTHING_SIGMA = 1.5;  % For Gaussian smoothing
        MAX_VALID_STRESS = 1e4;         % Pa, maximum reasonable stress value
        MIN_VALID_STRESS = -1e4;        % Pa, minimum reasonable stress value
        DEFAULT_TISSUE_THICKNESS = 10e-6;  % meters (10 µm)
    end
    
    methods (Static)
        function [stress_tensor, tension_map] = calculateStressAndTension(x, y, u, v, strain_tensor, E, nu, eta, dt, tissue_thickness)
            % Calculate stress tensor and tension from PIV and strain data
            % 
            % Inputs:
            %   x, y: Cell arrays of coordinate grids
            %   u, v: Cell arrays of displacements
            %   strain_tensor: Cell array of structs with xx, yy, xy components
            %   E: Young's modulus (Pa)
            %   nu: Poisson ratio
            %   eta: Viscosity (Pa·s)
            %   dt: Time step (seconds)
            %   tissue_thickness: Tissue thickness in meters (optional)
            %
            % Outputs:
            %   stress_tensor: Cell array of structs with stress components (Pa)
            %   tension_map: Cell array of tension values (N/m)
            
            % Set default tissue thickness if not provided
            if nargin < 10
                tissue_thickness = TissueStressAnalysis2.DEFAULT_TISSUE_THICKNESS;
            end
            
            % Validate inputs
            if ~iscell(strain_tensor)
                error('strain_tensor must be a cell array');
            end
            
            % Initialize outputs
            n_timepoints = length(strain_tensor);
            stress_tensor = cell(1, n_timepoints);
            tension_map = cell(1, n_timepoints);
            
            % Previous values for first iteration
            prev_exx = zeros(size(strain_tensor{1}.xx));
            prev_eyy = zeros(size(strain_tensor{1}.xx));
            prev_exy = zeros(size(strain_tensor{1}.xx));
            prev_sxx = zeros(size(strain_tensor{1}.xx));
            prev_syy = zeros(size(strain_tensor{1}.xx));
            prev_sxy = zeros(size(strain_tensor{1}.xx));
            
            % Process each timepoint
            for t = 1:n_timepoints
                % Get current strain components
                exx = strain_tensor{t}.xx;
                eyy = strain_tensor{t}.yy;
                exy = strain_tensor{t}.xy;
                
                % Clean data
                [exx_clean, eyy_clean, exy_clean] = TissueStressAnalysis2.cleanData(exx, eyy, exy);
                
                % Calculate stress components
                [sxx, syy, sxy] = TissueStressAnalysis2.calculateStressComponents(...
                    exx_clean, eyy_clean, exy_clean, ...    % Current strains
                    prev_exx, prev_eyy, prev_exy, ...       % Previous strains
                    prev_sxx, prev_syy, prev_sxy, ...       % Previous stresses
                    E, nu, eta, dt);                        % Material parameters
                
                % Store results
                stress_tensor{t} = struct('xx', sxx, 'yy', syy, 'xy', sxy);
                
                % Calculate surface tension (N/m)
                tension_map{t} = max((sxx + syy)/2, 0) * tissue_thickness;
                
                % Update previous values
                prev_exx = exx_clean;
                prev_eyy = eyy_clean;
                prev_exy = exy_clean;
                prev_sxx = sxx;
                prev_syy = syy;
                prev_sxy = sxy;
            end
        end
        
        function [sxx, syy, sxy] = calculateStressComponents(exx, eyy, exy, ...
                prev_exx, prev_eyy, prev_exy, ...
                prev_sxx, prev_syy, prev_sxy, ...
                E, nu, eta, dt)
            % Calculate strain rates
            dexx_dt = (exx - prev_exx) / dt;
            deyy_dt = (eyy - prev_eyy) / dt;
            dexy_dt = (exy - prev_exy) / dt;
            
            % Effective elastic modulus for plane stress
            E_eff = E / (1 - nu^2);
            
            % Maxwell model implementation
            % ds/dt + s/eta = E * de/dt
            
            % Calculate stress rates
            dsxx_dt = E_eff * (dexx_dt + nu*deyy_dt) - prev_sxx/eta;
            dsyy_dt = E_eff * (nu*dexx_dt + deyy_dt) - prev_syy/eta;
            dsxy_dt = E_eff * (1-nu) * dexy_dt - prev_sxy/eta;
            
            % Update stresses using forward Euler integration
            sxx = prev_sxx + dsxx_dt * dt;
            syy = prev_syy + dsyy_dt * dt;
            sxy = prev_sxy + dsxy_dt * dt;
            
            % Clip to reasonable values
            sxx = min(max(sxx, TissueStressAnalysis2.MIN_VALID_STRESS), ...
                TissueStressAnalysis2.MAX_VALID_STRESS);
            syy = min(max(syy, TissueStressAnalysis2.MIN_VALID_STRESS), ...
                TissueStressAnalysis2.MAX_VALID_STRESS);
            sxy = min(max(sxy, TissueStressAnalysis2.MIN_VALID_STRESS), ...
                TissueStressAnalysis2.MAX_VALID_STRESS);
        end
        
        function [exx_clean, eyy_clean, exy_clean] = cleanData(exx, eyy, exy)
            % Replace Inf and NaN with zeros
            exx_clean = exx;
            eyy_clean = eyy;
            exy_clean = exy;
            
            exx_clean(isinf(exx) | isnan(exx)) = 0;
            eyy_clean(isinf(eyy) | isnan(eyy)) = 0;
            exy_clean(isinf(exy) | isnan(exy)) = 0;
        end
        
        function visualizeStressMap(stress_tensor, tension_map, timepoint)
            figure('Position', [100 100 1200 400]);
            
            % Plot stress magnitude
            subplot(1,3,1);
            stress_mag = sqrt(stress_tensor{timepoint}.xx.^2 + ...
                stress_tensor{timepoint}.yy.^2 + ...
                stress_tensor{timepoint}.xy.^2);
            imagesc(stress_mag);
            colorbar;
            title('Stress Magnitude (Pa)');
            axis equal tight;
            
            % Plot tension
            subplot(1,3,2);
            imagesc(tension_map{timepoint});
            colorbar;
            title('Surface Tension (N/m)');
            axis equal tight;
            
            % Plot principal stress directions
            subplot(1,3,3);
            sxx = stress_tensor{timepoint}.xx;
            syy = stress_tensor{timepoint}.yy;
            sxy = stress_tensor{timepoint}.xy;
            
            [theta, lambda1, ~] = TissueStressAnalysis2.calculatePrincipalStresses(sxx, syy, sxy);
            
            % Downsample for visualization
            scale = 5;
            [X, Y] = meshgrid(1:scale:size(sxx,2), 1:scale:size(sxx,1));
            theta_down = theta(1:scale:end, 1:scale:end);
            lambda1_down = lambda1(1:scale:end, 1:scale:end);
            
            quiver(X, Y, cos(theta_down).*lambda1_down, sin(theta_down).*lambda1_down, 0.5);
            title('Principal Stress Directions');
            axis equal tight;
        end
        
        function [theta, lambda1, lambda2] = calculatePrincipalStresses(sxx, syy, sxy)
            % Calculate principal stresses and directions
            theta = 0.5 * atan2(2*sxy, sxx-syy);
            S = (sxx + syy)/2;
            R = sqrt(((sxx - syy)/2).^2 + sxy.^2);
            lambda1 = S + R;  % Maximum principal stress
            lambda2 = S - R;  % Minimum principal stress
        end
        
        function plotStressTimeSeries(stress_tensor, tension_map, dt)
            % Plot stress time series with proper time units
            n_timepoints = length(stress_tensor);
            time_vec = (0:(n_timepoints-1)) * dt/60;  % Convert to minutes
            
            % Initialize arrays
            mean_tension = zeros(1, n_timepoints);
            mean_sxx = zeros(1, n_timepoints);
            mean_syy = zeros(1, n_timepoints);
            mean_sxy = zeros(1, n_timepoints);
            mean_stress_mag = zeros(1, n_timepoints);
            
            std_tension = zeros(1, n_timepoints);
            std_sxx = zeros(1, n_timepoints);
            std_syy = zeros(1, n_timepoints);
            std_sxy = zeros(1, n_timepoints);
            std_stress_mag = zeros(1, n_timepoints);
            
            % Calculate statistics
            for t = 1:n_timepoints
                mean_tension(t) = mean(tension_map{t}(:), 'omitnan');
                std_tension(t) = std(tension_map{t}(:), 'omitnan');
                
                mean_sxx(t) = mean(stress_tensor{t}.xx(:), 'omitnan');
                mean_syy(t) = mean(stress_tensor{t}.yy(:), 'omitnan');
                mean_sxy(t) = mean(stress_tensor{t}.xy(:), 'omitnan');
                
                std_sxx(t) = std(stress_tensor{t}.xx(:), 'omitnan');
                std_syy(t) = std(stress_tensor{t}.yy(:), 'omitnan');
                std_sxy(t) = std(stress_tensor{t}.xy(:), 'omitnan');
                
                stress_mag = sqrt(stress_tensor{t}.xx.^2 + ...
                    stress_tensor{t}.yy.^2 + ...
                    stress_tensor{t}.xy.^2);
                mean_stress_mag(t) = mean(stress_mag(:), 'omitnan');
                std_stress_mag(t) = std(stress_mag(:), 'omitnan');
            end
            
            % Create figure
            figure('Position', [100 100 1200 800]);
            
            % Plot panels
            subplot(3,2,1);
            errorbar(time_vec, mean_tension, std_tension, 'b-', 'LineWidth', 2);
            xlabel('Time (minutes)');
            ylabel('Surface Tension (N/m)');
            title('Tissue Surface Tension');
            grid on;
            
            subplot(3,2,2);
            errorbar(time_vec, mean_stress_mag, std_stress_mag, 'k-', 'LineWidth', 2);
            xlabel('Time (minutes)');
            ylabel('Stress (Pa)');
            title('Total Stress Magnitude');
            grid on;
            
            subplot(3,2,3);
            errorbar(time_vec, mean_sxx, std_sxx, 'r-', 'LineWidth', 2);
            hold on;
            errorbar(time_vec, mean_syy, std_syy, 'g-', 'LineWidth', 2);
            xlabel('Time (minutes)');
            ylabel('Normal Stress (Pa)');
            title('Normal Stresses');
            legend('\sigma_{xx}', '\sigma_{yy}');
            grid on;
            
            subplot(3,2,4);
            errorbar(time_vec, mean_sxy, std_sxy, 'm-', 'LineWidth', 2);
            xlabel('Time (minutes)');
            ylabel('Shear Stress (Pa)');
            title('Shear Stress \tau_{xy}');
            grid on;
            
            subplot(3,2,5);
            d_tension = diff(mean_tension)./diff(time_vec);
            d_stress = diff(mean_stress_mag)./diff(time_vec);
            plot(time_vec(2:end), d_tension, 'b-', 'LineWidth', 2);
            hold on;
            plot(time_vec(2:end), d_stress, 'k--', 'LineWidth', 2);
            xlabel('Time (minutes)');
            ylabel('Rate (N/m/min, Pa/min)');
            title('Rate of Change');
            legend('Tension', 'Stress Magnitude');
            grid on;
            
            subplot(3,2,6);
            axis off;
            stats_text = {
                'Summary Statistics:';
                sprintf('Mean Tension: %.2f ± %.2f pN/um', mean(mean_tension)*1e6, mean(std_tension)*1e6);
                sprintf('Mean ?xx: %.2f ± %.2f Pa', mean(mean_sxx), mean(std_sxx));
                sprintf('Mean ?yy: %.2f ± %.2f Pa', mean(mean_syy), mean(std_syy));
                sprintf('Mean ?xy: %.2f ± %.2f Pa', mean(mean_sxy), mean(std_sxy));
                sprintf('Mean Total Stress: %.2f ± %.2f Pa', mean(mean_stress_mag), mean(std_stress_mag));
                };
            text(0.1, 0.8, stats_text, 'FontSize', 10, 'VerticalAlignment', 'top');
            
            sgtitle('Tissue Stress Analysis Over Time', 'FontSize', 14);
        end
        
        function measurements = calculateMechanicalMeasurements(stress_tensor, tension_map)
            n_timepoints = length(stress_tensor);
            
            % Initialize temporal arrays
            temporal = struct();
            temporal.tension = zeros(1, n_timepoints);
            temporal.pressure = zeros(1, n_timepoints);
            temporal.stress_magnitude = zeros(1, n_timepoints);
            temporal.normal_stress_xx = zeros(1, n_timepoints);
            temporal.normal_stress_yy = zeros(1, n_timepoints);
            temporal.shear_stress = zeros(1, n_timepoints);
            
            % Calculate temporal values
            for t = 1:n_timepoints
                temporal.tension(t) = mean(tension_map{t}(:), 'omitnan');
                
                sxx = stress_tensor{t}.xx(:);
                syy = stress_tensor{t}.yy(:);
                sxy = stress_tensor{t}.xy(:);
                
temporal.pressure(t) = -mean((sxx + syy)/2, 'omitnan');
                temporal.normal_stress_xx(t) = mean(sxx, 'omitnan');
                temporal.normal_stress_yy(t) = mean(syy, 'omitnan');
                temporal.shear_stress(t) = mean(abs(sxy), 'omitnan');
                
                % Calculate stress magnitude
                stress_mag = sqrt(sxx.^2 + syy.^2 + sxy.^2);
                temporal.stress_magnitude(t) = mean(stress_mag, 'omitnan');
            end
            
            % Calculate mean values
            mean_values = struct();
            fields = fieldnames(temporal);
            for i = 1:length(fields)
                field = fields{i};
                mean_values.(field) = mean(temporal.(field), 'omitnan');
            end
            
            % Calculate spatial maps
            [rows, cols] = size(stress_tensor{1}.xx);
            spatial = struct();
            spatial.tension_map = zeros(rows, cols);
            spatial.pressure_map = zeros(rows, cols);
            spatial.stress_magnitude_map = zeros(rows, cols);
            
            for t = 1:n_timepoints
                spatial.tension_map = spatial.tension_map + tension_map{t}/n_timepoints;
                
                % Calculate pressure map
                current_pressure = -(stress_tensor{t}.xx + stress_tensor{t}.yy)/2;
                spatial.pressure_map = spatial.pressure_map + current_pressure/n_timepoints;
                
                % Calculate stress magnitude map
                current_stress_mag = sqrt(stress_tensor{t}.xx.^2 + ...
                    stress_tensor{t}.yy.^2 + ...
                    stress_tensor{t}.xy.^2);
                spatial.stress_magnitude_map = spatial.stress_magnitude_map + ...
                    current_stress_mag/n_timepoints;
            end
            
            % Combine all measurements
            measurements = struct(...
                'mean_values', mean_values, ...
                'temporal', temporal, ...
                'spatial', spatial);
            
            % Display summary
            fprintf('\nMechanical Measurements Summary:\n');
            fprintf('--------------------------------\n');
            fprintf('Mean Tension: %.2f ± %.2f pN/um\n', ...
                mean_values.tension*1e6, std(temporal.tension, 'omitnan')*1e6);
            fprintf('Mean Pressure: %.2f ± %.2f Pa\n', ...
                mean_values.pressure, std(temporal.pressure, 'omitnan'));
            fprintf('Mean Stress Magnitude: %.2f ± %.2f Pa\n', ...
                mean_values.stress_magnitude, std(temporal.stress_magnitude, 'omitnan'));
            fprintf('Mean Normal Stress (xx): %.2f ± %.2f Pa\n', ...
                mean_values.normal_stress_xx, std(temporal.normal_stress_xx, 'omitnan'));
            fprintf('Mean Normal Stress (yy): %.2f ± %.2f Pa\n', ...
                mean_values.normal_stress_yy, std(temporal.normal_stress_yy, 'omitnan'));
            fprintf('Mean Shear Stress: %.2f ± %.2f Pa\n', ...
                mean_values.shear_stress, std(temporal.shear_stress, 'omitnan'));
        end
        
        function exportMeasurementsToCSV(measurements, output_folder)
            % Export measurements to CSV files
            % 
            % Inputs:
            %   measurements: struct from calculateMechanicalMeasurements
            %   output_folder: folder path to save CSV files
            
            if ~exist(output_folder, 'dir')
                mkdir(output_folder);
            end
            
            % 1. Export mean values
            mean_values_table = struct2table(measurements.mean_values);
            writetable(mean_values_table, fullfile(output_folder, 'mean_values.csv'));
            
            % 2. Export temporal data
            temporal_table = struct2table(measurements.temporal);
            writetable(temporal_table, fullfile(output_folder, 'temporal_evolution.csv'));
            
            % 3. Export spatial maps
            % Reshape 2D maps to tables
            writematrix(measurements.spatial.tension_map, ...
                fullfile(output_folder, 'spatial_tension.csv'));
            writematrix(measurements.spatial.pressure_map, ...
                fullfile(output_folder, 'spatial_pressure.csv'));
            writematrix(measurements.spatial.stress_magnitude_map, ...
                fullfile(output_folder, 'spatial_stress_magnitude.csv'));
            
            fprintf('Files exported to folder: %s\n', output_folder);
            fprintf('Generated files:\n');
            fprintf('- mean_values.csv\n');
            fprintf('- temporal_evolution.csv\n');
            fprintf('- spatial_tension.csv\n');
            fprintf('- spatial_pressure.csv\n');
            fprintf('- spatial_stress_magnitude.csv\n');
        end

function exportMeanValuesPerTimepointMasked(stress_tensor, tension_map, output_file, dt)
    % Export mean values per timepoint, excluding masked regions (zeros)
    
    if ~iscell(stress_tensor) || ~iscell(tension_map)
        error('stress_tensor and tension_map must be cell arrays');
    end
    
    if nargin < 4
        dt = 1;
    end
    
    n_timepoints = length(stress_tensor);
    
    mean_sxx = zeros(n_timepoints, 1);
    mean_syy = zeros(n_timepoints, 1);
    mean_sxy = zeros(n_timepoints, 1);
    mean_tension = zeros(n_timepoints, 1);
    mean_stress_mag = zeros(n_timepoints, 1);
    mean_pressure = zeros(n_timepoints, 1);
    mean_max_principal = zeros(n_timepoints, 1);
    mean_min_principal = zeros(n_timepoints, 1);
    
    time_vec = ((0:n_timepoints-1) * dt / 60)';  % Time in minutes
    
    for t = 1:n_timepoints
        sxx = stress_tensor{t}.xx;
        syy = stress_tensor{t}.yy;
        sxy = stress_tensor{t}.xy;
        tension = tension_map{t};
        
        % Create a mask that excludes zero-tension regions
        mask = tension ~= 0;
        
        % Compute derived quantities
        stress_mag = sqrt(sxx.^2 + syy.^2 + sxy.^2);
        pressure = -(sxx + syy)/2;
        [~, lambda1, lambda2] = TissueStressAnalysis.calculatePrincipalStresses(sxx, syy, sxy);
        
        % Only include nonzero tension pixels in averages
        mean_tension(t) = mean(tension(mask), 'omitnan');
        mean_sxx(t) = mean(sxx(mask), 'omitnan');
        mean_syy(t) = mean(syy(mask), 'omitnan');
        mean_sxy(t) = mean(sxy(mask), 'omitnan');
        mean_stress_mag(t) = mean(stress_mag(mask), 'omitnan');
        mean_pressure(t) = mean(pressure(mask), 'omitnan');
        mean_max_principal(t) = mean(lambda1(mask), 'omitnan');
        mean_min_principal(t) = mean(lambda2(mask), 'omitnan');
    end
    
    % Table output
    timepoint_data = table((0:n_timepoints-1)', time_vec, ...
        mean_sxx, mean_syy, mean_sxy, ...
        mean_tension*1e6, mean_stress_mag, ...
        mean_pressure, mean_max_principal, mean_min_principal, ...
        'VariableNames', {'Timepoint', 'Time_Minutes', ...
        'Mean_Stress_XX', 'Mean_Stress_YY', 'Mean_Stress_XY', ...
        'Mean_Tension', 'Mean_Stress_Magnitude', 'Mean_Pressure', ...
        'Mean_Max_Principal_Stress', 'Mean_Min_Principal_Stress'});
    
    % Write to file
    [output_dir, ~, ~] = fileparts(output_file);
    if ~isempty(output_dir) && ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    writetable(timepoint_data, output_file);
    
    fprintf('? Exported masked mean values to: %s\n', output_file);
end


    end
end