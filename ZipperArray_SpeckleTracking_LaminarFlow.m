%close all;
% Define simulation parameters
fx = 5e6;          % Center frequency: 5 MHz
c = 1540;          % Speed of sound in tissue (m/s)
wavelength = c / fx; 
gap = 0.01e-3;         
N_elements = 2;       

% ------------------------------
% index, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, material, width, height, dummy, dummy, dummy
% ------------------------------
% Using two rows (2 elements) data 

data = [... 
    1, -0.3e-3 - gap/2, (6e-3+gap)/2 - (2e-3) - gap, 0, -0.3e-3 - gap/2, (-6e-3-gap)/2, 0, -gap/2, (-6e-3-gap)/2, 0, -gap/2, (-6e-3-gap)/2 + (2e-3), 0, 1, 0.3e-3, 6e-3+gap, 0, 0, 0;...
    2, -0.3e-3 - gap/2, (6e-3+gap)/2, 0, -gap/2, (6e-3+gap)/2, 0, -gap/2, (-6e-3-gap)/2 + 2e-3 + gap, 0, -0.3e-3 - gap/2, (6e-3+gap)/2 - (2e-3), 0, 1, 0.3e-3, 6e-3+gap, 0, 0, 0];

% Calculate element centers (x, y, z) using four corners averaged
element_centers = zeros(N_elements, 3);
for i = 1:N_elements
    x_coords = [data(i,2), data(i,5), data(i,8), data(i,11)];
    y_coords = [data(i,3), data(i,6), data(i,9), data(i,12)];
    z_coords = [data(i,4), data(i,7), data(i,10), data(i,13)];
    element_centers(i,:) = [mean(x_coords), mean(y_coords), mean(z_coords)];
end
 
% Define odd and even element indices
odd_idx = 1:2:N_elements; 
even_idx = 2:2:N_elements;

% Print initial element centers using fixed formatting
fprintf('Initial odd element center(s):\n');
for i = odd_idx
    fprintf('  %.6f  %.6f  %.6f\n', element_centers(i,1), element_centers(i,2), element_centers(i,3));
end
fprintf('Initial even element center(s):\n');
for i = even_idx
    fprintf('  %.6f  %.6f  %.6f\n', element_centers(i,1), element_centers(i,2), element_centers(i,3));
end

% Calculate elevational pitch (only y)
odd_y = mean(element_centers(odd_idx,2)); 
even_y = mean(element_centers(even_idx,2));
y_distance = abs(odd_y - even_y) * 1e3; 
fprintf('Elevational pitch (y) = %.4f mm\n', y_distance);

% ------------------------------
% Grid definition (common for all iterations)
% ------------------------------
x_size = 10;        
z_size = 13;      
dx = 0.02;         
x = -x_size/2:dx:x_size/2;
z = 0:dx:z_size;
[X, Z] = meshgrid(x, z);

% ------------------------------
% Simulation parameters (common for all iterations)
% ------------------------------
num_scatterers = 3000;
y_size = 5;          % Scatterers distribution in y (the range)
sigma_y = 0.3;       % elevational beamwidth (mm)
p_comp = 1; 
     
% r = 0.9e-3;         % 血管徑向距離 [0.9,1.6,2.6,4.4] m
R = 5;                % 血管半徑 mm
Vmax = 10;            % 最大速度（mm/s）
dt = 1 / 100;         % 幀率 100 fps
% v_profile = Vmax * (1 - (r./R).^2);  % 單位 mm/s

num_frames = 50;            % Number of frames for dynamic simulation      
% y_move_speed = v_profile * dt * 1e-3;    % Movement per frame (m)
window_size = 128;          % Window size for correlation calculation
step_size = 4;              % Step size for correlation windows

% Create array to store max correlation coefficients for each frame across all iterations
num_iterations = 1;
all_max_cc = zeros(num_iterations, num_frames);
all_V_rep = zeros(1, num_iterations);

% Create figure for line plot
fig_line = figure('Name', 'Average Max CC vs Frame');

% Main loop for iterations
for iter = 1:num_iterations
    fprintf('\n=============================================\n');
    fprintf('Starting Iteration %d/%d\n', iter, num_iterations);
    fprintf('=============================================\n');
    
    % Generate new random scatterers for this iteration
    x_scatter = x_size * rand(num_scatterers, 1) - x_size/2;
    z_scatter = z_size * rand(num_scatterers, 1);
    y_scatter = y_size * rand(num_scatterers, 1) - y_size/2;
    amplitudes = randn(num_scatterers, 1) .* (1 + 0.5*rand(num_scatterers, 1));

    % vessel center
    vessel_cx = 0;
    vessel_cz = z_size / 2;
    % each scatter's radial position
    r = sqrt((x_scatter - vessel_cx).^2 + (z_scatter - vessel_cz).^2);
    % laminar flow 
    v_profile = Vmax * (1 - (r ./ R).^2);
    v_profile(r > R) = 0; % ensure only consider the scatterer in the vessel
    V_rep_iter = mean(v_profile(r<=R));    % 算平均速度 mm/s
    all_V_rep(iter) = V_rep_iter;
    % y direction displacement
    y_move_speed = v_profile * dt * 1e-3;
    
    % ------------------------------
    % (1) Reference image (Odd Tx/Odd Rx)
    % ------------------------------
    ref_centers = element_centers;
    composite_image_ref = zeros(size(X));
    % Use odd elements for transmit and receive
    for tx = odd_idx
        for rx = odd_idx
            sub_image = zeros(size(X));
            for s = 1:num_scatterers
                tx_delay = sqrt((x_scatter(s) - ref_centers(tx,1))^2 + ...
                                (y_scatter(s) - ref_centers(tx,2))^2 + ...
                                (z_scatter(s) - ref_centers(tx,3))^2) / c;
                rx_delay = sqrt((x_scatter(s) - ref_centers(rx,1))^2 + ...
                                (y_scatter(s) - ref_centers(rx,2))^2 + ...
                                (z_scatter(s) - ref_centers(rx,3))^2) / c;
                total_delay = tx_delay + rx_delay;
                
                r_tx = [x_scatter(s) - ref_centers(tx,1), ...
                        y_scatter(s) - ref_centers(tx,2), ...
                        z_scatter(s) - ref_centers(tx,3)];
                r_rx = [x_scatter(s) - ref_centers(rx,1), ...
                        y_scatter(s) - ref_centers(rx,2), ...
                        z_scatter(s) - ref_centers(rx,3)];
                directivity_tx = r_tx(3) / norm(r_tx);
                directivity_rx = r_rx(3) / norm(r_rx);
                directivity_orig = directivity_tx * directivity_rx;
                directivity = directivity_orig^(1 - p_comp); 
                
                elev_weight = exp(- (y_scatter(s)^2) / (2 * sigma_y^2) );
                
                sigma_x = 0.3;
                sigma_z = 0.3;
                % Fixed phase: removing random component
                phase = 2*pi*fx*total_delay;
                psf = amplitudes(s) * directivity * elev_weight * ...
                      exp(-(((X - x_scatter(s)).^2)/(2*sigma_x^2) + ((Z - z_scatter(s)).^2)/(2*sigma_z^2))) .* ...
                      cos(2*pi*fx/c * (Z - z_scatter(s)) + phase);
                sub_image = sub_image + psf;
            end
            composite_image_ref = composite_image_ref + sub_image;
        end
    end
    env_ref = abs(hilbert(composite_image_ref));
    env_norm_ref = env_ref / max(env_ref(:));
    ref_image = max(20*log10(env_norm_ref), -30);

    % ------------------------------
    % (2) Dynamic simulation: Even Elements (Even Tx/Even Rx)
    % ------------------------------
    
    % Copy initial element centers to moving_centers for simulation
    moving_centers = element_centers;
    
    % Display progress only in the first iteration
    if iter == 1
        % % Create a VideoWriter object to record the animation (only for first iteration)
        v = VideoWriter('even_animation.avi');
        v.FrameRate = 2;
        open(v);

        % Create a figure to display frames (only for first iteration)
        fig_frames = figure('Name','Dynamic Even-element Image (Even Tx/Even Rx)');
    end
    
    % Create cell array to store frame images (only for first iteration)
    if iter == 1
        frameImages = cell(1, num_frames);
    end
    
    for frame = 1:num_frames
        % Update scatter moving: move in negative y direction
        y_scatter = y_scatter + y_move_speed;
       
        % Print element centers only in the first iteration
        if iter == 1
            fprintf('Frame %d, current even element center(s):\n', frame);
            for i = even_idx
                fprintf('  %.6f  %.6f  %.6f\n', moving_centers(i,1), moving_centers(i,2), moving_centers(i,3));
            end
        end
        
        composite_image_even = zeros(size(X));
        % Use even elements for transmit and receive
        for tx = even_idx
            for rx = even_idx
                sub_image = zeros(size(X));
                for s = 1:num_scatterers
                    tx_delay = sqrt((x_scatter(s) - moving_centers(tx,1))^2 + ...
                                    (y_scatter(s) - moving_centers(tx,2))^2 + ...
                                    (z_scatter(s) - moving_centers(tx,3))^2) / c;
                    rx_delay = sqrt((x_scatter(s) - moving_centers(rx,1))^2 + ...
                                    (y_scatter(s) - moving_centers(rx,2))^2 + ...
                                    (z_scatter(s) - moving_centers(rx,3))^2) / c;
                    total_delay = tx_delay + rx_delay;
                    
                    r_tx = [x_scatter(s) - moving_centers(tx,1), ...
                            y_scatter(s) - moving_centers(tx,2), ...
                            z_scatter(s) - moving_centers(tx,3)];
                    r_rx = [x_scatter(s) - moving_centers(rx,1), ...
                            y_scatter(s) - moving_centers(rx,2), ...
                            z_scatter(s) - moving_centers(rx,3)];
                    directivity_tx = r_tx(3) / norm(r_tx);
                    directivity_rx = r_rx(3) / norm(r_rx);
                    directivity_orig = directivity_tx * directivity_rx;
                    directivity = directivity_orig^(1 - p_comp);
                    
                    elev_weight = exp(- (y_scatter(s)^2) / (2 * sigma_y^2) );
                    
                    sigma_x = 0.3;
                    sigma_z = 0.3;
                    % Fixed phase: removing random component
                    phase = 2*pi*fx*total_delay;
                    psf = amplitudes(s) * directivity * elev_weight * ...
                          exp(-(((X - x_scatter(s)).^2)/(2*sigma_x^2) + ((Z - z_scatter(s)).^2)/(2*sigma_z^2))) .* ...
                          cos(2*pi*fx/c * (Z - z_scatter(s)) + phase);
                    sub_image = sub_image + psf;
                end
                composite_image_even = composite_image_even + sub_image;
            end
        end
        
        env_even = abs(hilbert(composite_image_even));
        env_norm_even = env_even / max(env_even(:));
        even_image = max(20*log10(env_norm_even), -30);
        
        % Calculate correlation map and use maximum correlation coefficient as indicator
        corr_map = calculateCorrelationMap(ref_image, even_image, window_size, step_size);
        max_cc = max(corr_map(:));
        %mean_cc = mean(corr_map(:));
        
        % Store max correlation coefficient for this frame and iteration
        all_max_cc(iter, frame) = max_cc;
        %all_max_cc(iter, frame) = mean_cc;
        
        
        % Print max CC for each frame only in the first iteration
        if iter == 1
            fprintf('Frame %d: Max Correlation Coefficient = %.4f\n', frame, max_cc);
            %fprintf('Frame %d: Mean Correlation Coefficient = %.4f\n', frame, mean_cc);
        end
        
        % Display and record frames only for the first iteration
        if iter == 1
            figure(fig_frames);
            imagesc(x, z, even_image);
            colormap(gray); colorbar;
            xlabel('Lateral Distance (mm)'); ylabel('Axial Distance (mm)');
            title(sprintf('Frame %d: Even Elements, Max CC = %.4f', frame, max_cc));
            %title(sprintf('Frame %d: Even Elements, Max CC = %.4f', frame, mean_cc));
            drawnow;
            
            % Capture the frame and write to video file
            frame_img = getframe(gcf);
            writeVideo(v, frame_img);
            
            % Save the current frame image in cell array
            frameImages{frame} = even_image;
        end
    end
    % figure each scatterer velocity
    if iter == 1
        % 3D
        figure('Name', '3D Scatterer Velocity Profile');
        x_plot_mm = (x_scatter - vessel_cx) ; 
        z_plot_mm = (z_scatter - vessel_cz);
        scatter3(x_plot_mm, z_plot_mm, v_profile, 30, v_profile, 'filled');
        xlabel('相對於血管中心的 X 位置 (mm)');
        ylabel('相對於血管中心的 Z 位置 (mm)');
        zlabel('速度 (mm/s)');
        title(sprintf('三維層流速度剖面 (Vmax = %.1f mm/s, R = %.1f mm)', Vmax, R));
        cb = colorbar; % 顯示顏色條
        ylabel(cb, '速度 (mm/s)');
        axis equal; % 使得 x-z 平面上的圓形看起來更像圓形
        view(-35, 45); % 調整視角以便觀察，您可以嘗試不同的角度
        grid on;
    end
    % Print iteration progress
    fprintf('Completed iteration %d/%d\n', iter, num_iterations);
end

% Close video writer (only created in the first iteration)
if exist('v', 'var') && isvalid(v)
    close(v);
end

% ------------------------------
% Calculate and display average max CC for each frame
% ------------------------------
avg_max_cc = mean(all_max_cc, 1);
fprintf('\n=============================================\n');
fprintf('Average Max Correlation Coefficient for Each Frame (across %d iterations):\n', num_iterations);
for frame = 1:num_frames
    fprintf('Frame %d: Avg Max CC = %.4f\n', frame, avg_max_cc(frame));
end
% ------------------------------
% Find max CC and frame --> using mean velocity
% ------------------------------
V_rep = mean(all_V_rep);    % mm/s
[peak_cc_value, peak_frame_index] = max(avg_max_cc);
fprintf('平均最大相關係數在第 %d 幀達到峰值，其值為: %.4f\n', peak_frame_index, peak_cc_value);
displacement_at_peak_cc_mm = Vmax * dt * peak_frame_index;
fprintf('在相關係數峰值幀 (第 %d 幀)，平均流速為 %.3f mm/s，累積位移為: %.4f mm\n', peak_frame_index, V_rep, displacement_at_peak_cc_mm);

% ------------------------------
% Plot average max CC vs frame number with accumulated displacement
% ------------------------------
figure(fig_line);
plot(1:num_frames, avg_max_cc, 'o-', 'MarkerSize', 3, 'MarkerFaceColor', 'b', 'LineWidth', 2);
xlabel('Frame Number');
ylabel('Average Max Correlation Coefficient');
title(sprintf('Average Max CC vs Frame (%d Iterations)', num_iterations));
grid on;

% Add secondary x-axis for accumulated displacement

accumulated_displacement = (1:num_frames) .* y_move_speed * 1e3; % Convert to mm
ax1 = gca;
% ax2 = axes('Position', get(ax1, 'Position'), 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none', 'XColor', 'r', 'YColor', 'none');
% line(accumulated_displacement, avg_max_cc, 'Parent', ax2, 'Color', 'r', 'LineStyle', '-', 'Marker', 'x', 'LineWidth', 2);
% xlabel(ax2, 'Accumulated Displacement (mm)', 'Color', 'r');
% xlim(ax2, [min(accumulated_displacement), max(accumulated_displacement)]);
xlim(ax1, [1, num_frames]);

% Create a legend
legend(ax1, 'Avg Max CC vs Frame');
% legend(ax2, 'Avg Max CC vs Displacement');

% ------------------------------
% Plot individual frame images from the first iteration
% ------------------------------
% if exist('frameImages', 'var')
%     for k = 1:num_frames
%         figure('Name', sprintf('Frame %d', k));
%         imagesc(x, z, frameImages{k});
%         colormap(gray); colorbar;
%         title(sprintf('Frame %d', k));
%         xlabel('Lateral (mm)');
%         ylabel('Axial (mm)');
%         xticks([]); yticks([]);
%     end
% end

% ------------------------------
% Plot CC values for all iterations and average
% ------------------------------
figure('Name', 'CC Comparison Across Iterations');
hold on;
for iter = 1:num_iterations
    plot(1:num_frames, all_max_cc(iter, :), 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
end
plot(1:num_frames, avg_max_cc, 'o-', 'MarkerSize', 3, 'MarkerFaceColor', 'b', 'LineWidth', 2, 'Color', 'r');
xlabel('Frame Number');
ylabel('Max Correlation Coefficient');
title('Max CC for All Iterations with Average');
grid on;
legend('Individual Iterations', 'Average', 'Location', 'southwest');

% ------------------------------
% compare with odd image and maxCC even image
% ------------------------------
best_even = frameImages{peak_frame_index};  
ref_img  = ref_image;

figure('Name','Reference vs. Best-Matched Frame');
subplot(1,2,1);
imagesc(x, z, ref_img);
colormap(gray); axis image off;
title('Odd-element Reference');
subplot(1,2,2);
imagesc(x, z, best_even);
colormap(gray); axis image off;
title(sprintf('Even-element Frame %d (Max CC)', peak_frame_index));

diff_img = best_even - ref_img;
figure('Name','Difference Image');
imagesc(x, z, diff_img);
colormap(gray); axis image off;
title('Even – Odd Difference');
colorbar;
%% cc function
function [corr_map, best_pos] = calculateCorrelationMap(orig_img, moved_img, window_size, step_size)
    [h, w] = size(orig_img);
    num_windows_y = floor((h - window_size)/step_size) + 1;
    num_windows_x = floor((w - window_size)/step_size) + 1;
    corr_map = zeros(num_windows_y, num_windows_x);
    % init 
    max_cc = -Inf;
    best_pos = [1, 1];

    for i = 1:num_windows_y
        y_start = (i-1)*step_size + 1;
        y_end = y_start + window_size - 1;
        for j = 1:num_windows_x
            x_start = (j-1)*step_size + 1;
            x_end = x_start + window_size - 1;
            win1 = orig_img(y_start:y_end, x_start:x_end);
            win2 = moved_img(y_start:y_end, x_start:x_end);
            R = corrcoef(win1(:), win2(:));
            corr_map(i,j) = R(1,2);
            if R(1,2) > max_cc
                max_cc = R(1,2);
                best_pos = [x_start, y_start];  % matched block 左上角座標
            end
        end
    end
end
