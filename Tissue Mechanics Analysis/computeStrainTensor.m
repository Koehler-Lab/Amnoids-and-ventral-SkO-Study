function [strain_tensors] = computeStrainTensor(x, y, u, v)
    % Initialize strain tensor storage
    n_frames = numel(x);  % Assuming x is a cell array
    strain_tensors = cell(1, n_frames);
    
    % Compute strain tensor for each frame
    for frame = 1:n_frames
        % Extract frame data
        current_x = x{frame};
        current_y = y{frame};
        current_u = u{frame};
        current_v = v{frame};
        
        % Create masks for valid data
        valid_u = ~isnan(current_u) & ~isinf(current_u);
        valid_v = ~isnan(current_v) & ~isinf(current_v);
        
        % Initialize strain tensor components
        strain_tensor = struct(...
            'xx', zeros(size(current_u)), ...   % Normal strain in x-direction
            'yy', zeros(size(current_v)), ...   % Normal strain in y-direction
            'xy', zeros(size(current_u)) ...    % Shear strain
        );
        
        % Compute spatial derivatives using central difference method
        [ny, nx] = size(current_x);
        
        % Strain computation with boundary handling
        for i = 2:(ny-1)
            for j = 2:(nx-1)
                % Skip computation if any required points contain NaN or Inf
                if ~all(valid_u(i, [j-1,j,j+1])) || ~all(valid_u([i-1,i,i+1], j)) || ...
                   ~all(valid_v(i, [j-1,j,j+1])) || ~all(valid_v([i-1,i,i+1], j))
                    strain_tensor.xx(i,j) = NaN;
                    strain_tensor.yy(i,j) = NaN;
                    strain_tensor.xy(i,j) = NaN;
                    continue;
                end
                
                % Spatial steps
                dx = current_x(i,j+1) - current_x(i,j-1);
                dy = current_y(i+1,j) - current_y(i-1,j);
                
                % Skip if spatial steps are too small
                if abs(dx) < eps || abs(dy) < eps
                    strain_tensor.xx(i,j) = NaN;
                    strain_tensor.yy(i,j) = NaN;
                    strain_tensor.xy(i,j) = NaN;
                    continue;
                end
                
                % Displacement gradients
                % Normal strain xx: du/dx
                strain_tensor.xx(i,j) = (current_u(i,j+1) - current_u(i,j-1)) / (2*dx);
                
                % Normal strain yy: dv/dy
                strain_tensor.yy(i,j) = (current_v(i+1,j) - current_v(i-1,j)) / (2*dy);
                
                % Shear strain xy: 0.5 * (du/dy + dv/dx)
                strain_tensor.xy(i,j) = 0.5 * (...
                    (current_u(i+1,j) - current_u(i-1,j)) / (2*dy) + ...
                    (current_v(i,j+1) - current_v(i,j-1)) / (2*dx) ...
                );
            end
        end
        
        % Handle edge boundaries
        strain_tensor = handleBoundaries(strain_tensor, current_x, current_y, current_u, current_v);
        
        % Apply smoothing
        strain_tensor = smoothStrainTensor(strain_tensor);
        
        % Store strain tensor
        strain_tensors{frame} = strain_tensor;
    end
end

function [strain_tensor] = handleBoundaries(strain_tensor, x, y, u, v)
    % Handle boundary conditions
    [ny, nx] = size(x);
    
    % Create masks for valid data at boundaries
    valid_u = ~isnan(u) & ~isinf(u);
    valid_v = ~isnan(v) & ~isinf(v);
    
    % Top and bottom edges (forward/backward differences)
    for j = 1:nx
        % Top edge
        if all(valid_u(1:2,j)) && all(valid_v(1:2,j))
            strain_tensor.xx(1,j) = (u(1,j) - u(2,j)) / (x(2,j) - x(1,j));
            strain_tensor.yy(1,j) = (v(1,j) - v(2,j)) / (y(2,j) - y(1,j));
        else
            strain_tensor.xx(1,j) = NaN;
            strain_tensor.yy(1,j) = NaN;
        end
        
        % Bottom edge
        if all(valid_u(end-1:end,j)) && all(valid_v(end-1:end,j))
            strain_tensor.xx(end,j) = (u(end,j) - u(end-1,j)) / (x(end,j) - x(end-1,j));
            strain_tensor.yy(end,j) = (v(end,j) - v(end-1,j)) / (y(end,j) - y(end-1,j));
        else
            strain_tensor.xx(end,j) = NaN;
            strain_tensor.yy(end,j) = NaN;
        end
    end
    
    % Left and right edges
    for i = 1:ny
        % Left edge
        if all(valid_u(i,1:2)) && all(valid_v(i,1:2))
            strain_tensor.xx(i,1) = (u(i,1) - u(i,2)) / (x(i,2) - x(i,1));
            strain_tensor.yy(i,1) = (v(i,1) - v(i,2)) / (y(i,2) - y(i,1));
        else
            strain_tensor.xx(i,1) = NaN;
            strain_tensor.yy(i,1) = NaN;
        end
        
        % Right edge
        if all(valid_u(i,end-1:end)) && all(valid_v(i,end-1:end))
            strain_tensor.xx(i,end) = (u(i,end) - u(i,end-1)) / (x(i,end) - x(i,end-1));
            strain_tensor.yy(i,end) = (v(i,end) - v(i,end-1)) / (y(i,end) - y(i,end-1));
        else
            strain_tensor.xx(i,end) = NaN;
            strain_tensor.yy(i,end) = NaN;
        end
    end
    
    % Shear strain at boundaries
    strain_tensor.xy([1,end],:) = NaN;
    strain_tensor.xy(:,[1,end]) = NaN;
end

function [strain_tensor] = smoothStrainTensor(strain_tensor)
    % Apply Gaussian smoothing to reduce noise, handling NaN values
    strain_tensor.xx = nanfilt(strain_tensor.xx, 1);
    strain_tensor.yy = nanfilt(strain_tensor.yy, 1);
    strain_tensor.xy = nanfilt(strain_tensor.xy, 1);
end

function smoothed = nanfilt(A, sigma)
    % Custom filtering function that handles NaN values
    B = A;
    nanMask = isnan(A);
    
    % Only proceed if there are some NaN values but not all
    if any(nanMask(:)) && ~all(nanMask(:))
        % Replace NaNs with mean of non-NaN neighbors
        B(nanMask) = 0;
        
        % Apply Gaussian filter
        smoothed = imgaussfilt(B, sigma, 'FilterSize', 5);
        
        % Restore NaN values where they were originally
        smoothed(nanMask) = NaN;
    else
        % If no NaNs or all NaNs, just apply regular filtering
        smoothed = imgaussfilt(A, sigma, 'FilterSize', 5);
    end
end

% % Cell array input (common in PIV data)
% strain_tensors = computeMultiFrameStrainTensor(...
%     x, y, u_original, v_original, ...
%     'time_step', 0.1, ...  % Time between frames
%     'smoothing', true, ...
%     'interpolation', true ...
% );
% 
% % 3D matrix input
% strain_tensors = computeMultiFrameStrainTensor(...
%     x_3d, y_3d, u_3d, v_3d ...
% );