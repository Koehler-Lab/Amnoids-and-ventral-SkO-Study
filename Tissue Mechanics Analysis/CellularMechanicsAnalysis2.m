classdef CellularMechanicsAnalysis2
    properties (Constant)
        % Physical and numerical constants
        BOLTZMANN_CONSTANT = 1.380649e-23;  % J/K
        ROOM_TEMPERATURE = 293.15;          % K
    end
    
    methods (Static)
        
        % Utility: Truncated Normal Distribution Sampling
        function [val] = truncatedNormal(mu, sigma, min_val, max_val)
            % Generate truncated normal random variable
            while true
                val = normrnd(mu, sigma);
                if val >= min_val && val <= max_val
                    break;
                end
            end
        end
        
        % Main analysis function
function [results] = analyzeCellularMechanics(strain_tensor, segmentation_mask, varargin)
    fprintf('Starting cellular mechanics analysis...\n');
    
    % Parse input arguments
    fprintf('Parsing input arguments...\n');
    p = inputParser;
    addParameter(p, 'base_priors', [], @isstruct);
    parse(p, varargin{:});
    opts = p.Results;

    % Extract cell metrics from segmentation
    fprintf('Analyzing cell segmentation...\n');
    cell_metrics = CellularMechanicsAnalysis2.analyzeCellSegmentation(segmentation_mask);
    fprintf('Cell segmentation complete - Found %d cells\n', cell_metrics.num_cells);

    % Define base priors if not provided
    fprintf('Defining priors...\n');
    if isempty(opts.base_priors)
        base_priors = CellularMechanicsAnalysis2.defineDefaultPriors();
    else
        base_priors = opts.base_priors;
    end
    fprintf('Priors defined\n');

    % Add debug prints right before and after the function call
    fprintf('Debug: About to call modifyBayesianEstimationWithCellMetrics\n');
    fprintf('Debug: strain_tensor class: %s\n', class(strain_tensor));
    fprintf('Debug: cell_metrics class: %s\n', class(cell_metrics));
    fprintf('Debug: base_priors class: %s\n', class(base_priors));
    
    try
        % Modify Bayesian estimation with cellular metrics
        fprintf('Starting Bayesian estimation...\n');
        results = CellularMechanicsAnalysis2.modifyBayesianEstimationWithCellMetrics(...
            strain_tensor, cell_metrics, base_priors);
        fprintf('Debug: Bayesian estimation completed\n');
    catch ME
        fprintf('Error occurred:\n');
        fprintf('Message: %s\n', ME.message);
        fprintf('Stack:\n');
        disp(ME.stack);
        rethrow(ME);
    end
    
    fprintf('Analysis complete\n');
end
        % Cell Segmentation Analysis
function [metrics] = analyzeCellSegmentation(segmentation_mask)
    % Get boundaries
    boundaries = segmentation_mask == 255;
    
    % Dilate boundaries slightly to ensure closure
    se = strel('disk', 1);
    boundaries_closed = imdilate(boundaries, se);
    
    % Fill in cells
    filled = imfill(~boundaries_closed, 'holes');
    
    % Erode to compensate for dilation
    filled = imerode(filled, se);
    
    % Remove small objects and border objects
    filled_clean = bwareaopen(filled, 100);  % Remove small objects
    filled_no_border = imclearborder(filled_clean);
    
    % Label individual cells
    [labeled_cells, num_cells] = bwlabel(filled_no_border, 8);  % Using 8-connectivity
    
    % Visualize steps for debugging
    figure('Name', 'Segmentation Steps');
    
    subplot(2,3,1);
    imshow(segmentation_mask);
    title('Original Mask');
    
    subplot(2,3,2);
    imshow(boundaries_closed);
    title('Closed Boundaries');
    
    subplot(2,3,3);
    imshow(filled);
    title('After Filling');
    
    subplot(2,3,4);
    imshow(filled_clean);
    title('After Cleaning');
    
    subplot(2,3,5);
    imshow(filled_no_border);
    title('After Border Removal');
    
    subplot(2,3,6);
    imshow(label2rgb(labeled_cells));
    title(['Labeled Cells: ' num2str(num_cells)]);
    
    fprintf('Found %d cells\n', num_cells);
    
    if num_cells == 0
        error('No cells detected. Please check segmentation parameters.');
    end
    
    % Get cell properties
    cell_props = regionprops(labeled_cells, 'Area', 'Perimeter', 'Centroid', ...
        'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Eccentricity');
    
    % Extract and calculate all metrics
    cell_areas = [cell_props.Area];
    cell_perimeters = [cell_props.Perimeter];
    
    metrics = struct();
    metrics.total_area = sum(cell_areas);
    metrics.mean_cell_area = mean(cell_areas);
    metrics.median_cell_area = median(cell_areas);
    metrics.cell_area_std = std(cell_areas);
    metrics.area_variation_coefficient = metrics.cell_area_std / metrics.mean_cell_area;
    
    metrics.roundness = 4 * pi * cell_areas ./ (cell_perimeters.^2);
    metrics.mean_roundness = mean(metrics.roundness);
    
    major_axes = [cell_props.MajorAxisLength];
    minor_axes = [cell_props.MinorAxisLength];
    metrics.aspect_ratios = major_axes ./ minor_axes;
    metrics.mean_aspect_ratio = mean(metrics.aspect_ratios);
    
    cell_orientations = [cell_props.Orientation];
    metrics.orientation_distribution = histcounts(cell_orientations, 18);
    metrics.mean_orientation = mean(cell_orientations);
    
    metrics.neighbor_distances = CellularMechanicsAnalysis2.computeNeighborDistances(cell_props);
    metrics.num_cells = num_cells;
    metrics.cell_density = num_cells / metrics.total_area;
    
    metrics.eccentricity = [cell_props.Eccentricity];
    metrics.mean_eccentricity = mean(metrics.eccentricity);
end

            function [neighbor_distances] = computeNeighborDistances(cell_props)
                % Calculate distances between cell centroids
                centroids = cat(1, cell_props.Centroid);
                num_cells = size(centroids, 1);

                % Initialize distance matrix
                neighbor_distances = zeros(num_cells, num_cells);

                % Compute pairwise distances
                for i = 1:num_cells
                    for j = i+1:num_cells
                        dist = norm(centroids(i,:) - centroids(j,:));
                        neighbor_distances(i,j) = dist;
                        neighbor_distances(j,i) = dist;
                    end
                end
            end

        function [mechanical_properties, property_distributions] = estimateMechanicalPropertiesWithCellContext(mcmc_samples, cell_metrics)
        % Estimate mechanical properties with cellular context

        % Flatten and reshape samples
        samples_flat = reshape(mcmc_samples, [], 3);

        % Compute basic statistical properties
        mechanical_properties = struct();
        param_names = {'E', 'nu', 'eta'};

        for i = 1:3
            param = samples_flat(:, i);
            mechanical_properties.(param_names{i}) = struct(...
                'mean', mean(param), ...
                'median', median(param), ...
                'std', std(param), ...
                'ci', prctile(param, [2.5, 97.5]), ...
                'samples', param ...
            );
        end
        
        

        % Compute property distributions
        property_distributions = struct();
        for i = 1:3
            [counts, edges] = histcounts(samples_flat(:,i), 50, 'Normalization', 'probability');
            property_distributions.(param_names{i}) = struct(...
                'counts', counts, ...
                'edges', edges ...
            );
        end

        % Add cellular context interpretation
        mechanical_properties.cellular_context = struct(...
            'heterogeneity_score', cell_metrics.area_variation_coefficient, ...
            'shape_factor', mean(cell_metrics.aspect_ratios), ...
            'packing_density', cell_metrics.cell_density ...
        );
        end

        function visualizeResults(results)
        % Create comprehensive visualization of analysis results
        figure('Name', 'Cellular Mechanics Analysis', 'Position', [100, 100, 1500, 1000]);

        % 1. Mechanical Property Distributions
        param_names = {'E', 'nu', 'eta'};
        param_labels = {'Young''s Modulus (Pa)', 'Poisson Ratio', 'Viscosity (Pa·s)'};

        for i = 1:3
            subplot(2,3,i);
            histogram(results.mechanical_properties.(param_names{i}).samples, 50, ...
                'Normalization', 'probability', 'FaceColor', 'b', 'EdgeColor', 'k');
            title(sprintf('%s Distribution', param_labels{i}));
            xlabel(param_labels{i});
            ylabel('Probability');
            grid on;
        end

        % 2. Prior Modification Impact
        subplot(2,3,4);
        modification_info = results.prior_modification;
        modification_factors = [
            modification_info.area_variation_impact, ...
            modification_info.shape_factor, ...
            modification_info.packing_factor
        ];
        bar(modification_factors);
        title('Prior Modification Factors');
        xticks(1:3);
        xticklabels({'Area Variation', 'Shape Factor', 'Packing Factor'});
        ylabel('Factor Value');
        grid on;

        % 3. MCMC Diagnostics
        subplot(2,3,5);
        bar(results.sampling_diagnostics.r_hat);
        title('MCMC Convergence (R-hat)');
        xticks(1:3);
        xticklabels(param_names);
        ylabel('R-hat Value');
        grid on;

        % 4. Cellular Context Summary
        subplot(2,3,6);
        context = results.mechanical_properties.cellular_context;
        metrics = [
            context.heterogeneity_score, ...
            context.shape_factor, ...
            context.packing_density
        ];
        bar(metrics);
        title('Cellular Mechanical Context');
        xticks(1:3);
        xticklabels({'Heterogeneity', 'Shape Factor', 'Density'});
        ylabel('Metric Value');
        grid on;

        % Adjust layout
        sgtitle('Cellular Mechanics Analysis Results', 'FontSize', 14);

        % Try to save figure
        try
            saveas(gcf, 'cellular_mechanics_analysis.png');
        catch
            warning('Could not save visualization');
        end
        end

function residuals = computeResiduals(params, strain_tensor)
    % Extract parameters
    E = params(1);    % Young's modulus
    nu = params(2);   % Poisson ratio
    eta = params(3);  % Viscosity
    
    % Clean strain tensor data
    xx = strain_tensor.xx;
    yy = strain_tensor.yy;
    xy = strain_tensor.xy;
    
    % Replace Inf and NaN with finite values or remove them from calculation
    max_valid_strain = 1.0;  % Maximum reasonable strain value
    
    % Function to clean data
    function clean_data = cleanStrainData(data)
        clean_data = data;
        % Replace Inf with max valid value
        clean_data(isinf(data)) = sign(data(isinf(data))) * max_valid_strain;
        % Replace NaN with 0 or nearby valid values
        clean_data(isnan(data)) = 0;
        % Clip extreme values
        clean_data = min(max(clean_data, -max_valid_strain), max_valid_strain);
    end
    
    % Clean all components
    xx_clean = cleanStrainData(xx);
    yy_clean = cleanStrainData(yy);
    xy_clean = cleanStrainData(xy);
    
    % Calculate with cleaned data
    % Elastic contribution
    elastic_xx = xx_clean / E;
    elastic_yy = yy_clean / E;
    elastic_xy = xy_clean / (2*(1 + nu));
    
    % Viscous contribution if time data available
    if isfield(strain_tensor, 'time')
        t = strain_tensor.time;
        viscous_xx = xx_clean * t / eta;
        viscous_yy = yy_clean * t / eta;
        viscous_xy = xy_clean * t / eta;
    else
        viscous_xx = zeros(size(xx_clean));
        viscous_yy = zeros(size(yy_clean));
        viscous_xy = zeros(size(xy_clean));
    end
    
    % Total theoretical strains
    theoretical_xx = elastic_xx + viscous_xx;
    theoretical_yy = elastic_yy + viscous_yy;
    theoretical_xy = elastic_xy + viscous_xy;
    
    % Compute residuals only for valid data points
    valid_points = ~isnan(xx) & ~isinf(xx) & ...
                  ~isnan(yy) & ~isinf(yy) & ...
                  ~isnan(xy) & ~isinf(xy);
    
    residuals_xx = (xx_clean - theoretical_xx) .* valid_points;
    residuals_yy = (yy_clean - theoretical_yy) .* valid_points;
    residuals_xy = (xy_clean - theoretical_xy) .* valid_points;
    
    % Combine residuals
    residuals = [residuals_xx(:); residuals_yy(:); residuals_xy(:)];
    
    % Only use valid points in final residuals
    residuals = residuals(isfinite(residuals));
end

    % Define Default Priors
        function [base_priors, modification_info] = defineDefaultPriors(varargin)
            % Define default mechanical property priors with optional modifications

            % Default prior specifications
            base_priors.E = struct(...
                'mean', 1000, ...    % Young's Modulus (Pa)
                'std', 300, ...     % Standard deviation
                'min', 100, ...     % Minimum value
                'max', 5000 ...     % Maximum value
            );

            base_priors.nu = struct(...
                'mean', 0.48, ...   % Poisson Ratio
                'std', 0.002, ...    % Standard deviation
                'min', 0.46, ...    % Minimum value
                'max', 0.499 ...    % Maximum value
            );

            base_priors.eta = struct(...
                'mean', 500, ...  % Viscosity (Pa·s)
                'std', 150, ...   % Standard deviation
                'min', 10, ...   % Minimum value
                'max', 5000 ...   % Maximum value
            );

            % Optional modification of priors
            if nargin > 0 && isstruct(varargin{1})
                modification = varargin{1};

                % Modify priors if specific fields are provided
                if isfield(modification, 'E')
                    if isfield(modification.E, 'mean')
                        base_priors.E.mean = modification.E.mean;
                    end
                    if isfield(modification.E, 'std')
                        base_priors.E.std = modification.E.std;
                    end
                end

                % Similar modifications for nu and eta
                if isfield(modification, 'nu')
                    if isfield(modification.nu, 'mean')
                        base_priors.nu.mean = modification.nu.mean;
                    end
                    if isfield(modification.nu, 'std')
                        base_priors.nu.std = modification.nu.std;
                    end
                end

                if isfield(modification, 'eta')
                    if isfield(modification.eta, 'mean')
                        base_priors.eta.mean = modification.eta.mean;
                    end
                    if isfield(modification.eta, 'std')
                        base_priors.eta.std = modification.eta.std;
                    end
                end
            end

            % Prepare modification information
            modification_info = struct(...
                'E_mean', base_priors.E.mean, ...
                'E_std', base_priors.E.std, ...
                'nu_mean', base_priors.nu.mean, ...
                'nu_std', base_priors.nu.std, ...
                'eta_mean', base_priors.eta.mean, ...
                'eta_std', base_priors.eta.std ...
            );
        end

        % Modify Bayesian Estimation with Cellular Metrics
function [modified_estimation] = modifyBayesianEstimationWithCellMetrics(strain_tensor, cell_metrics, base_priors)
    fprintf('Debug 1: Entered modifyBayesianEstimationWithCellMetrics\n');
    
    % Check inputs
    fprintf('Debug 2: Checking inputs\n');
    fprintf('Strain tensor fields: %s\n', strjoin(fieldnames(strain_tensor), ', '));
    fprintf('Cell metrics fields: %s\n', strjoin(fieldnames(cell_metrics), ', '));
    fprintf('Base priors fields: %s\n', strjoin(fieldnames(base_priors), ', '));
    
    fprintf('Debug 3: About to create cellular informed priors\n');
    try
        [refined_priors, prior_modification_info] = CellularMechanicsAnalysis2.createCellularInformedPriors(...
            base_priors, cell_metrics);
        fprintf('Debug 4: Created cellular informed priors\n');
    catch ME
        fprintf('Error in createCellularInformedPriors: %s\n', ME.message);
        rethrow(ME);
    end
    
    fprintf('Debug 5: About to setup likelihood function\n');
    try
        likelihood_function = @(params) CellularMechanicsAnalysis2.computeCellularWeightedLikelihood(...
            params, strain_tensor, cell_metrics);
        fprintf('Debug 6: Created likelihood function\n');
    catch ME
        fprintf('Error in creating likelihood function: %s\n', ME.message);
        rethrow(ME);
    end
    
    fprintf('Debug 7: Starting MCMC sampling\n');
    [mcmc_samples, sampling_diagnostics] = CellularMechanicsAnalysis2.runCellularInformedMCMC(...
        refined_priors, likelihood_function, strain_tensor);
    
    fprintf('Debug 8: Estimating mechanical properties\n');
    [mechanical_properties, property_distributions] = CellularMechanicsAnalysis2.estimateMechanicalPropertiesWithCellContext(...
        mcmc_samples, cell_metrics);
    
    modified_estimation = struct(...
        'mechanical_properties', mechanical_properties, ...
        'property_distributions', property_distributions, ...
        'prior_modification', prior_modification_info, ...
        'sampling_diagnostics', sampling_diagnostics);
end

        % Create Cellular-Informed Priors
function [refined_priors, modification_info] = createCellularInformedPriors(base_priors, cell_metrics)
    fprintf('Starting prior refinement...\n');
    
    % Print input data
    fprintf('Base priors structure:\n');
    disp(base_priors);
    
    fprintf('Cell metrics structure:\n');
    disp(cell_metrics);
    
    try
        % Young's Modulus Prior Modification
        fprintf('Modifying E prior...\n');
        E_variation = cell_metrics.area_variation_coefficient;
        refined_priors.E = struct(...
            'mean', base_priors.E.mean * (1 + E_variation), ...
            'std', base_priors.E.std * (1 + E_variation), ...
            'min', base_priors.E.min, ...
            'max', base_priors.E.max ...
        );
        
        % Poisson Ratio Prior Modification
        fprintf('Modifying nu prior...\n');
        cell_shape_factor = 1 + cell_metrics.mean_eccentricity;
        refined_priors.nu = struct(...
            'mean', base_priors.nu.mean * cell_shape_factor, ...
            'std', base_priors.nu.std * cell_shape_factor, ...
            'min', base_priors.nu.min, ...
            'max', base_priors.nu.max ...
        );
        
        % Viscosity Prior Modification
        fprintf('Modifying eta prior...\n');
        cell_packing_factor = min(1 / cell_metrics.cell_density,10);
        refined_priors.eta = struct(...
            'mean', base_priors.eta.mean * min(cell_packing_factor, 10), ...  % Cap the scaling
            'std', base_priors.eta.std * min(cell_packing_factor, 10), ...
            'min', base_priors.eta.min, ...
            'max', min(base_priors.eta.max * cell_packing_factor, 5000) ... % Cap the maximum
        );
        
        % Store modification information
        fprintf('Creating modification info...\n');
        modification_info = struct(...
            'area_variation_impact', E_variation, ...
            'shape_factor', cell_shape_factor, ...
            'packing_factor', cell_packing_factor ...
        );
        
        fprintf('Prior refinement complete\n');
        
    catch ME
        fprintf('Error in createCellularInformedPriors: %s\n', ME.message);
        fprintf('Error occurred in: %s\n', ME.stack(1).name);
        rethrow(ME);
    end
end
        
        % Compute Cellular-Weighted Likelihood
        function [likelihood] = computeCellularWeightedLikelihood_v1(...
                params, strain_tensor, cell_metrics)
            % Compute likelihood with cellular metric weighting
            
            % Standard mechanical likelihood (placeholder - implement your specific likelihood)
            base_likelihood = CellularMechanicsAnalysis2.computeStandardLikelihood(...
                params, strain_tensor);
            
            % Cellular metric penalties
            cellular_penalties = CellularMechanicsAnalysis2.computeCellularPenalties(...
                cell_metrics);
            
            % Combined likelihood
            likelihood = base_likelihood + cellular_penalties;
            
            return;
        end
        
function [likelihood] = computeCellularWeightedLikelihood(params, strain_tensor, cell_metrics)
    % Compute likelihood with cellular metric weighting
    
    % Get base mechanical likelihood
    base_likelihood = CellularMechanicsAnalysis2.computeStandardLikelihood(...
        params, strain_tensor);
    
    % Get cellular penalties
    cellular_penalties = CellularMechanicsAnalysis2.computeCellularPenalties(...
        cell_metrics);
    
    % Combined likelihood with weight of 0.1 for penalties
    likelihood = base_likelihood + 0.1 * cellular_penalties;
end        
        % Compute Standard Likelihood (Placeholder)
function [log_likelihood] = computeStandardLikelihood(params, strain_tensor)
    
    % Very large sigma for more acceptance
    sigma = 1000.0;
    
    try
        % Get residuals and print their values
        residuals = CellularMechanicsAnalysis2.computeResiduals(params, strain_tensor);

        % Check for infinities or NaNs in residuals
        if any(isinf(residuals(:))) || any(isnan(residuals(:)))
            fprintf('Warning: Found Inf or NaN in residuals!\n');
            log_likelihood = -1e6;
            return;
        end
        
        % More numerically stable computation
        scaled_residuals = residuals(:) ./ sigma;
        sum_sq = sum(scaled_residuals.^2);
                
        % Compute log likelihood with protection against overflow
        log_likelihood = -0.5 * sum_sq;
        
        % Add prior terms
        nu_prior = -0.5 * ((params(2) - 0.48) / 0.01)^2;
        E_prior = -0.5 * ((params(1) - 1500) / 300)^2;
        eta_prior = -0.5 * ((params(3) - 1500) / 300)^2;
        
        log_likelihood = log_likelihood + nu_prior + E_prior + eta_prior;
        
    catch ME
        fprintf('Error in likelihood calculation: %s\n', ME.message);
        fprintf('Error occurred in: %s, line %d\n', ME.stack(1).name, ME.stack(1).line);
        log_likelihood = -1e6;
    end
end

% Compute Cellular Penalties
        function [cellular_penalties] = computeCellularPenalties(cell_metrics)
            % Compute penalty terms based on cellular metrics
            
            % Area variation penalty
            area_variation_penalty = -log(1 + cell_metrics.area_variation_coefficient);
            
            % Shape irregularity penalty
            shape_penalty = -log(cell_metrics.mean_roundness);
            
            % Eccentricity influence
            eccentricity_penalty = -cell_metrics.mean_eccentricity;
            
            % Combined penalties
            cellular_penalties = area_variation_penalty + ...
                                 shape_penalty + ...
                                 eccentricity_penalty;
            
            return;
        end
        
        % Run Cellular-Informed MCMC
function [mcmc_samples, sampling_diagnostics] = runCellularInformedMCMC(refined_priors, likelihood_function, strain_tensor)
    % MCMC parameters
    n_chains = 1;  % Reduced to 1 chain for now
    samples_per_chain = 10000;
    burnin = 5000;
    
    % Preallocate arrays with correct dimensions
    effective_samples = samples_per_chain - burnin;
    mcmc_samples = zeros(n_chains, effective_samples, 3);  % Changed dimension
    acceptance_rates = zeros(n_chains, 1);
    
    % Run parallel chains
    for i = 1:n_chains
        fprintf('Starting chain %d/%d\n', i, n_chains);
        [chain_samples, acceptance_rate] = CellularMechanicsAnalysis2.runSingleAdaptiveMCMCChain(...
            refined_priors, likelihood_function, samples_per_chain, burnin);
        
        mcmc_samples(i,:,:) = chain_samples;  % This should now match dimensions
        acceptance_rates(i) = acceptance_rate;
    end
    
    % Compute convergence diagnostics
    r_hat = CellularMechanicsAnalysis2.computeConvergenceDiagnostics(mcmc_samples);
    
    % Return diagnostics
    sampling_diagnostics = struct(...
        'acceptance_rates', acceptance_rates, ...
        'mean_acceptance_rate', mean(acceptance_rates), ...
        'r_hat', r_hat ...
    );
end

function [chain_samples, acceptance_rate] = runSingleAdaptiveMCMCChain(refined_priors, likelihood_function, samples_per_chain, burnin)
    % Initialize
    current_params = [1500, 0.48, 1500];  % [E, nu, eta]
    chain_samples = zeros(samples_per_chain - burnin, 3);
    
    % Different scales for different parameters
    param_scales = zeros(1,3);
    param_scales(1) = 500;     % E: steps of 100 Pa
    param_scales(2) = 0.002;   % nu: tiny steps around 0.48
    param_scales(3) = 500;     % eta: steps of 100 Pa·s
    
    n_accepted = 0;
    
    for i = 1:samples_per_chain
        % Generate proposals separately for each parameter
        proposal = zeros(1,3);
        proposal(1) = current_params(1) + randn * param_scales(1);  % E
        proposal(2) = current_params(2) + randn * param_scales(2);  % nu
        proposal(3) = current_params(3) + randn * param_scales(3);  % eta
        
       % Print progress every 1000 iterations
        if mod(i, 10000) == 0
            fprintf('Iteration %d/%d (%.1f%%) - Accepted: %d (%.1f%%)\n', ...
                i, samples_per_chain, 100*i/samples_per_chain, ...
                n_accepted, 100*n_accepted/i);
        end

        
        % Check bounds separately
        in_bounds = true;
        if proposal(1) < refined_priors.E.min || proposal(1) > refined_priors.E.max || ...
           proposal(2) < refined_priors.nu.min || proposal(2) > refined_priors.nu.max || ...
           proposal(3) < refined_priors.eta.min || proposal(3) > refined_priors.eta.max
            in_bounds = false;
        end
        
        if in_bounds
            try
                current_ll = likelihood_function(current_params);
                proposal_ll = likelihood_function(proposal);
                
                if mod(i, 1000) == 0
                    fprintf('Current LL: %.4e, Proposal LL: %.4e\n', current_ll, proposal_ll);
                end
                
                % Metropolis acceptance
                log_ratio = proposal_ll - current_ll;
                if log(rand) < log_ratio
                    current_params = proposal;
                    n_accepted = n_accepted + 1;
                    if mod(i, 100) == 0
                        fprintf('Accepted!\n');
                    end
                end
                
            catch ME
                fprintf('Error in iteration %d: %s\n', i, ME.message);
                continue;
            end
        end
        
        if i > burnin
            chain_samples(i-burnin,:) = current_params;
        end
    end
    
    % Print final statistics
    acceptance_rate = n_accepted / samples_per_chain;
    fprintf('\nMCMC Final Statistics:\n');
    fprintf('Acceptance rate: %.2f%%\n', 100*acceptance_rate);
    fprintf('Final parameters: E=%.2f, nu=%.4f, eta=%.2f\n', ...
        current_params(1), current_params(2), current_params(3));
    fprintf('Parameter ranges:\n');
    fprintf('E: min=%.2f, max=%.2f, std=%.2f\n', ...
        min(chain_samples(:,1)), max(chain_samples(:,1)), std(chain_samples(:,1)));
    fprintf('eta: min=%.2f, max=%.2f, std=%.2f\n', ...
        min(chain_samples(:,3)), max(chain_samples(:,3)), std(chain_samples(:,3)));
end

        % Compute Convergence Diagnostics
    function [r_hat] = computeConvergenceDiagnostics(mcmc_samples)
        % Compute Gelman-Rubin convergence diagnostic
        [n_chains, n_samples, n_params] = size(mcmc_samples);

        r_hat = zeros(n_params, 1);

        for p = 1:n_params
            % Extract parameter samples
            param_samples = squeeze(mcmc_samples(:,:,p));

            % Calculate chain means and overall mean
            chain_means = mean(param_samples, 2);
            overall_mean = mean(chain_means);

            % Between-chain variance
            B = n_samples * var(chain_means);

            % Within-chain variance
            W = mean(var(param_samples, 0, 2));

            % Pooled variance estimate
            V_hat = ((n_samples - 1)/n_samples) * W + ((n_chains + 1)/(n_chains * n_samples)) * B;

            % Calculate R-hat
            r_hat(p) = sqrt(V_hat/W);
        end
    end

    end
end

% % Example Usage Function
% function [results] = runCellularMechanicsAnalysis2(strain_tensor, segmentation_mask)
%     % Wrapper function for easy analysis
%     
%     % Perform comprehensive cellular mechanics analysis
%     results = CellularMechanicsAnalysis2.analyzeCellularMechanics(...
%         strain_tensor, segmentation_mask);
%     
%     % Print key results summary
%     fprintf('\n=== Cellular Mechanics Analysis Results ===\n');
%     fprintf('Young''s Modulus: %.2f ± %.2f Pa\n', ...
%         results.mechanical_properties.E.mean, ...
%         results.mechanical_properties.E.std);
%     fprintf('Poisson Ratio: %.4f ± %.4f\n', ...
%         results.mechanical_properties.nu.mean, ...
%         results.mechanical_properties.nu.std);
%     fprintf('Viscosity: %.2f ± %.2f Pa·s\n', ...
%         results.mechanical_properties.eta.mean, ...
%         results.mechanical_properties.eta.std);
%     
%     % Print cellular context description
%     fprintf('\nTissue Behavior Description:\n%s\n', ...
%         results.mechanical_properties.cellular_context.tissue_behavior_description);
% end

% 
% % Assuming you have:
% % - strain_tensor: Computed strain tensor from displacement data
% % - segmentation_mask: Binary mask of cell boundaries
% 
% % Run the analysis
% results = runCellularMechanicsAnalysis2(strain_tensor, segmentation_mask);