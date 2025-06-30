% close all;
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
% Define vessel center: ROI (Region of Interest)
% ------------------------------
roi_size = 2.5;           % ROI size (mm)
half_size = roi_size / 2;
radius_tolerance = 0.1;   % vessel center (mm)
vessel_cx = 0;            % vessel center X coordinate
vessel_cz = z_size / 2;   % vessel center Z coordinate
x_min = vessel_cx - half_size;  % roi square
x_max = vessel_cx + half_size;
z_min = vessel_cz - half_size;
z_max = vessel_cz + half_size;
center_mask = (X >= x_min) & (X <= x_max) & (Z >= z_min) & (Z <= z_max);
if ~any(center_mask(:))
    warning('Vessel center ROI is empty. Check grid definition.');
end
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
dt = 1 / 100;         % 幀率 50/100/150/200 fps

num_frames = 50;            % Number of frames for dynamic simulation      
% y_move_speed = v_profile * dt * 1e-3;    % Movement per frame (m)
window_size = 128;          % Window size for correlation calculation
step_size = 4;              % Step size for correlation windows

% Create array to store max correlation coefficients for each frame across all iterations
num_iterations = 2;
all_center_cc = zeros(num_iterations, num_frames); % save vessel center cc 
all_V_rep = zeros(1, num_iterations);
measured_velocities_per_iteration = zeros(1, num_iterations);

% Create figure for line plot
% fig_line = figure('Name', 'Average Max CC vs Frame');

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

    % each scatter's radial position
    r = sqrt((x_scatter - vessel_cx).^2 + (z_scatter - vessel_cz).^2);
    % laminar flow 
    v_profile = Vmax * (1 - (r ./ R).^2);
    v_profile(r > R) = 0; % ensure only consider the scatterer in the vessel
    
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
        
        % find vessel center in ref_img & even_img 
        ref_pixels = ref_image(center_mask);
        even_pixels = even_image(center_mask);
        if ~isempty(ref_pixels) && ~isempty(even_pixels)
            C = corrcoef(ref_pixels, even_pixels);
            center_cc = C(1, 2);
        else
            center_cc = NaN; % 如果區域為空，則結果為 NaN (Not a Number)
        end
        all_center_cc(iter, frame) = center_cc;
        
        % Print max CC for each frame only in the first iteration
        if iter == 1
            fprintf('Frame %d: Center ROI Correlation = %.4f\n', frame, center_cc);
            %fprintf('Frame %d: Mean Correlation Coefficient = %.4f\n', frame, mean_cc);
        end
        
        % Display and record frames only for the first iteration
        if iter == 1
            figure(fig_frames);
            imagesc(x, z, even_image);
            colormap(gray); colorbar;
            xlabel('Lateral Distance (mm)'); ylabel('Axial Distance (mm)');
            title(sprintf('Frame %d: Even Elements, Center CC = %.4f', frame, center_cc));
            % 確認血管半徑跟 ROI 位置
            % hold on;
            % vessel_position = [vessel_cx - R, vessel_cz - R, 2*R, 2*R];
            % rectangle('Position', vessel_position, 'Curvature', [1,1], ...
            %           'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '-');
            % roi_position = [x_min, z_min, x_max - x_min, z_max - z_min];
            % rectangle('Position', roi_position, 'EdgeColor', 'y', 'LineWidth', 1.5, 'LineStyle', '--');
            % 
            % hold off;
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
%  vessel center's mean cc
%  ------------------------------
% 'omitnan' 會在計算平均值時，自動忽略任何 NaN 的數值
avg_center_cc = mean(all_center_cc, 1, 'omitnan');

fprintf('\n=============================================\n');
fprintf('Average Correlation Coefficient at the Vessel Center:\n');
for frame = 1:num_frames
    fprintf('Frame %d: Avg CC = %.4f\n', frame, avg_center_cc(frame));
end

% ------------------------------
%  find the vessel center's MaxCC and frame
%  ------------------------------
% 忽略第一幀 (因位移為0，理論上CC最高)，以找出由位移產生的峰值
if num_frames > 1
    [initial_peak_cc, initial_peak_frame] = max(avg_center_cc(2:end));
    initial_peak_frame = initial_peak_frame + 1;
else
    [initial_peak_cc, initial_peak_frame] = max(avg_center_cc);
end
% 二次曲線內插法，獲得更精確的峰值 
if initial_peak_frame > 1 && initial_peak_frame < num_frames
    % 取出初始峰值點和其左右相鄰點的 y 值
    y1 = avg_center_cc(initial_peak_frame - 1);
    y2 = avg_center_cc(initial_peak_frame);    % 這就是 initial_peak_cc
    y3 = avg_center_cc(initial_peak_frame + 1);
    
    % 拋物線擬合公式，計算相對於中間點(peak_frame)的偏移量 delta
    delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
    
    % 計算出帶有小數的、更精確的峰值幀數
    precise_peak_frame = initial_peak_frame + delta;
    
    % 我們也可以估算出這條拋物線真正的頂點高度，作為更精確的 peak_cc
    precise_peak_cc = y2 - 0.25 * (y1 - y3) * delta;
else
    % 如果峰值在邊界，則不使用內插，沿用整數結果
    precise_peak_frame = initial_peak_frame;
    precise_peak_cc = initial_peak_cc;
end
% 流速計算(mm/s)
peak_displacement = y_distance; 
time_to_peak = precise_peak_frame * dt;
measured_velocity = y_distance / time_to_peak;

fprintf('\n---------------------------------------------\n');
fprintf('Final Analysis for the Vessel Center:\n');
fprintf('---------------------------------------------\n');
fprintf('  -> Peak Correlation:         %.4f\n', precise_peak_cc);
fprintf('  -> Occurred at Frame:        %.4f\n', precise_peak_frame);
fprintf('  -> Time to Peak:             %.4f s\n', time_to_peak);
fprintf('  -> Implied Displacement:     %.4f mm (Equal to Elev. Pitch)\n', peak_displacement);
fprintf('  -> Measured Velocity:        %.4f mm/s\n', measured_velocity);
fprintf('---------------------------------------------\n');
fprintf('\n--- Verification ---\n');
fprintf('  Simulated Velocity at Center (Vmax): %.2f mm/s\n', Vmax);
fprintf('  Measured Velocity from CC Peak:    %.2f mm/s\n', measured_velocity);
fprintf('  Error Value:    %.2f \n', (Vmax-measured_velocity)/Vmax*100);

fprintf('--------------------\n');

% ------------------------------
%  mean cc v.s frame
%  ------------------------------
figure('Name', 'Correlation vs. Frame for Vessel Center');
plot(1:num_frames, avg_center_cc, 'o-', 'LineWidth', 2);
grid on;
title(sprintf('Average Correlation at Vessel Center (%d Iterations)', num_iterations));
xlabel('Frame Number');
ylabel('Average Correlation Coefficient');
legend('Vessel Center');

% ------------------------------
%  繪製所有模擬的CC曲線與平均值，檢查穩定性
%  ------------------------------
figure('Name', 'CC Comparison Across Iterations for Vessel Center');
hold on;
% 畫出每一次獨立模擬的曲線 (灰色)
for iter = 1:num_iterations
    plot(1:num_frames, all_center_cc(iter, :), 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
end
% 畫出平均曲線 (紅色粗體)
plot(1:num_frames, avg_center_cc, 'o-', 'MarkerSize', 3, 'LineWidth', 2, 'Color', 'r');
hold on;
plot(precise_peak_frame, precise_peak_cc, 'r*', 'MarkerSize', 10, 'LineWidth', 2); % 用一個紅星標出峰值
hold off;
result_text = sprintf('Measured Velocity: %.2f mm/s', measured_velocity);
annotation('textbox', [0.65, 0.8, 0.1, 0.1], 'String', result_text, 'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black');
grid on;
xlabel('Frame Number');
ylabel('Correlation Coefficient');
title('All Iteration CCs and Average at Vessel Center');
legend('Individual Iterations', 'Average', 'Location', 'southwest');

%  -----------------------------
%  save data
%  -----------------------------
filename = sprintf('Vmax_%.1f_FrameRate_%.2f_ROI_%.1f_iter_%d_results.mat', Vmax, dt, roi_size, num_iterations);
save(filename, 'avg_center_cc', 'all_center_cc', 'precise_peak_cc', 'precise_peak_frame', 'measured_velocity', 'Vmax', 'dt', 'y_distance', 'measured_velocities_per_iteration', 'mean_measured_velocity', 'std_measured_velocity');

fprintf('Results saved to file: %s\n', filename);
