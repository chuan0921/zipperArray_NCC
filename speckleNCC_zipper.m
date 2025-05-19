function speckleNCC_zipper()
    %% 初始化參數
    settings = struct(...
        'num_iterations', 15, ...        % 幀數
        'fx', 5e6, ...
        'c', 1540, ...
        'dt', 1/100, ...                 % 幀間時間 (秒)
        'Vmax', 25, ...                  % 最大速度 mm/s
        'num_scatterers', 5000, ...
        'vessel_radius', 5, ...
        'vessel_length', 14, ...
        'min_distance', 0.02, ...
        'vessel_r', 0.9, ...             % kernel 中心半徑
        'z_size', 13, ...                % 成像深度 13*10mm
        'x_size', 10, ...
        'dz', 0.05, ...
        'dx', 0.05, ...
        'zipper_num', 4, ...            % (先不用)Zipper 元件數
        'zipper_pitch', 3, ...          % mm
        'zipper_width', 0.3, ...
        'zipper_height', 6, ...
        'zipper_gap', 0.01 ...
    );
    % zipper_data = [...
    %     1, settings.zipper_gap/2, (settings.zipper_width + settings.zipper_gap)/2 - 2, 0, ...
    %        settings.zipper_gap/2, (-settings.zipper_height - settings.zipper_gap)/2, 0, ...
    %        settings.zipper_gap/2 + settings.zipper_width, (-settings.zipper_height - settings.zipper_gap)/2, 0, ...
    %        settings.zipper_width + settings.zipper_gap/2, (-settings.zipper_width - settings.zipper_gap)/2 + 2, 0, ...
    %        1, settings.zipper_width, settings.zipper_height + settings.zipper_gap, ...
    %        0, 0, 0;  % 共18個數據 + 1個 id
    %     2, settings.zipper_gap/2, (settings.zipper_width + settings.zipper_gap)/2, 0, ...
    %        settings.zipper_gap/2, (settings.zipper_width + settings.zipper_gap)/2 - 2, 0, ...
    %        settings.zipper_width + settings.zipper_gap/2, (-settings.zipper_width - settings.zipper_gap)/2 + 2, 0, ...
    %        settings.zipper_gap/2 + settings.zipper_width, (settings.zipper_width + settings.zipper_gap)/2, 0, ...
    %        1, settings.zipper_width, settings.zipper_width + settings.zipper_gap, ...
    %        0, 0, 0;
    % ];
    zipper_data = [ ...
    1, 0,    0.0,    0, ...     % P1_L (ele_center for Left)
       -0.005,  1.0,    0, ...  % P2_L 
       -0.005, -1.995,  0, ...  % P3_L 
       -2.995, -1.995,  0, ...  % P4_L 
       -2.995,  1.995,  0, ...  % P5_L 
        0,      0,      0; ...  % P6_L 

    2,  0,    3.01,   0, ...    % P1_R (ele_center for Right)
        0.005,  5.005,  0, ...  % P2_R 
        0.005,  2.01,   0, ...  % P3_R 
        2.995,  1.015,  0, ...  % P4_R 
        2.995,  5.005,  0, ...  % P5_R 
        0,      0,      0       % P6_R 
];
    rng(42);   % 固定scatter
    num_repeat = 1;    % 跑20次不同 speckle樣本
    [X, Z, x, z] = buildMeshgrid(settings);
    
    ncc_all = zeros(num_repeat, settings.num_iterations);
    displacement_all = zeros(num_repeat, settings.num_iterations);  % 只記錄 dy 位移
    velocity_all = zeros(num_repeat, settings.num_iterations);
    dy_peak = zeros(num_repeat, 1);
    
    for trial = 1:num_repeat
        fprintf('第 %d 次模擬...\n', trial);

        [sx, sy, sz, amp] = genScatters(settings);
        sx_iter = sx; sy_iter = sy; sz_iter = sz;

        num_sub_elements = size(zipper_data, 1);
        tx_list_left = 1:2:num_sub_elements;
        rx_list_left = 1:2:num_sub_elements;
        tx_list_right = 2:2:num_sub_elements;
        rx_list_right = 2:2:num_sub_elements;

        [original_image, moved_images] = generateAllImages(sx_iter, sy_iter, sz_iter, amp, tx_list_left, rx_list_left, tx_list_right, rx_list_right, settings, X, Z,zipper_data);
        [image_h, image_w] = size(original_image);
        % 自動調整kernel size
        frac_x = 0.6;  % 20% image size 
        frac_z = 0.6;  
        kernel_mm_x = frac_x * settings.x_size;  
        kernel_mm_z = frac_z * settings.z_size;  
    
        % 转成像素并 clamp
        kernel_w = min(round(kernel_mm_x / settings.dx), image_w);
        kernel_h = min(round(kernel_mm_z / settings.dz), image_h);
    
        % 居中放置
        cx = floor((image_w - kernel_w)/2) + 1;
        cz = floor((image_h - kernel_h)/2) + 1;
    
        % clamp 起点
        kx = max(1, min(cx, image_w - kernel_w + 1));
        kz = max(1, min(cz, image_h - kernel_h + 1));
        kernel_block = [kx, kz, kx + kernel_w - 1, kz + kernel_h - 1];
        matched_positions = cell(1, settings.num_iterations);
        %--- 計算 origin_image 跟每張 moved_images 的 NCC，並找出最高的那張
        for l = 1:settings.num_iterations
            [disp_vec, cc_max, matched_pos] = calculateMatchedBlock_3D(original_image, moved_images{l}, kernel_block, kernel_w, kernel_h, settings, l, x, z);
            if l == 1 % moved_images{1} 是 t=0 的影像，Y位移為0
                current_dy = 0;
            else % moved_images{l} for l>1, 粒子移動了 l-1 步
                current_dy = (l-1) * settings.Vmax * settings.dt;
            end
            displacement_all(trial, l) = disp_vec(2);
            velocity_all(trial, l) = disp_vec(2) / ((l-1) * settings.dt+ eps);
            ncc_all(trial, l) = cc_max;
            matched_positions{l} = matched_pos;
        end
        % 取最佳相關係數對應的影像
        [best_ncc, best_idx] = max(ncc_all(trial, :));
        best_frame = best_idx;
        best_dy    = displacement_all(trial, best_idx);
        best_vel   = velocity_all(trial, best_idx);
        fprintf('最佳匹配：第 %d 幀 (NCC=%.3f)，位移=%.3f mm，速度=%.3f mm/s\n', ...
                best_frame, best_ncc, best_dy, best_vel);
        dy_peak(trial) = best_dy;

        disp(displacement_all);
        disp(ncc_all);
        disp("左元件座標：");
        disp(zipper_data(1, 2:4));
        disp("右元件座標：");
        disp(zipper_data(2, 2:4));
        fprintf("真實位移 per frame: %.2f mm\n", settings.Vmax * settings.dt);
        fprintf("幀數: %d\n", settings.num_iterations);
        %% 驗證scatter流動位置
        fprintf('\n--- 散射粒子 Y 座標驗證 (Trial %d, Best Frame %d) ---\n', trial, best_frame);
        [~, sy_at_best_idx_step, ~, ~] = updateScatterers(sx, sy, sz, amp, settings, best_idx-1);
        y_center_right_element = zipper_data(2,3); % 右元件的 Y 中心
        
        fprintf('在 Best Frame (對應 %d 個 dt 的移動後):\n', best_idx-1);
        fprintf('  右元件 Y 中心目標: %.4f mm\n', y_center_right_element);
        fprintf('  此時散射粒子 Y 座標平均值: %.4f mm\n', mean(sy_at_best_idx_step));
        fprintf('  此時散射粒子 Y 座標標準差: %.4f mm\n', std(sy_at_best_idx_step));
        fprintf('  此時散射粒子 Y 座標最小值: %.4f mm\n', min(sy_at_best_idx_step));
        fprintf('  此時散射粒子 Y 座標最大值: %.4f mm\n', max(sy_at_best_idx_step));
        fprintf('  此時散射粒子數量: %d\n', numel(sy_at_best_idx_step));
    end
    %% NCC vs y displacement
    figure;
    hold on;
    for trial = 1:num_repeat
        plot(displacement_all(trial,:), ncc_all(trial,:), 'Color', [0.7 0.9 1]);
    end
    plot(mean(displacement_all, 1, 'omitnan'), mean(ncc_all, 1, 'omitnan'), '-b', 'LineWidth', 2);
    xlabel('Y Displacement (mm)'); ylabel('NCC'); grid on;
    title('NCC vs Y Displacement'); hold off;
    %% Speckle 移動影片（新增 kernel 與 matched block 顯示）
    v = VideoWriter('speckle_flow.avi');
    v.FrameRate = 5;
    open(v);
    fig = figure('Name', 'Speckle Flow Frame');
    for i = 1:settings.num_iterations
        imagesc(x, z, moved_images{i}); colormap gray; axis image; colorbar;
        title(sprintf('Frame %d', i)); hold on;
        rectangle('Position', [(x(kx)), (z(kz)), kernel_w * settings.dx, kernel_h * settings.dz], 'EdgeColor', 'g', 'LineWidth', 1.5);
        if ~isempty(matched_positions{i})
            rect_x = matched_positions{i}(1);
            rect_z = matched_positions{i}(2);
            rectangle('Position', [rect_x, rect_z, kernel_w * settings.dx, kernel_h * settings.dz], ...
                      'EdgeColor', 'm', 'LineStyle', '--', 'LineWidth', 1.2);
        end
        hold off;
        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end
    close(v);
    close(fig);
    implay('speckle_flow.avi');
    
     %% 原始影像與 kernel 標示
    figure('Name', 'Original Image');
    imagesc(x, z, original_image); colormap gray; colorbar;
    xlabel('Lateral (mm)'); ylabel('Axial (mm)'); axis image;
    rectangle('Position', [(x(kx)), (z(kz)), kernel_w * settings.dx, kernel_h * settings.dz], 'EdgeColor', 'r', 'LineWidth', 1.5);
end
function [X, Z, x, z] = buildMeshgrid(settings)
    x = -settings.x_size/2 : settings.dx : settings.x_size/2;
    z = 0 : settings.dz : settings.z_size;
    [X, Z] = meshgrid(x, z);
end

function [sx, sy, sz, amp] = genScatters(s)
    max_trials = s.num_scatterers * 10;  % 嘗試上限
    points = zeros(max_trials, 3);
    count = 0;

    while count < s.num_scatterers
        candidate = [(rand()-0.5)*s.x_size, ...
                     (rand()-0.5)*s.zipper_height, ...
                     rand()*s.z_size];

        if count == 0 || all(vecnorm(points(1:count,:) - candidate, 2, 2) > s.min_distance)
            count = count + 1;
            points(count, :) = candidate;
        end
    end

    sx = points(1:count, 1);
    sy = points(1:count, 2);
    sz = points(1:count, 3);
    amp = sqrt(-2 * log(rand(count, 1)));  % Rayleigh 分布
end

function [sx, sy, sz, amp] = updateScatterers(sx, sy, sz, amp, s ,step)
    % 固定 scatter 模式 + y 方向流動
    if nargin < 6
        step = 1;  % 預設每次只移動一幀
    end
    % vy = s.Vmax * ones(size(sy));          % 每個點都以最大速度往 y 方向移動之後要改層流
    vy = s.Vmax;  
    sy = sy + vy * s.dt * step;            % 更新位置：y = y + v * dt * step
    
    % 保留在成像區域內的 scatterers
    inside_idx = (abs(sx) <= s.x_size/2) & (abs(sy) <= s.zipper_height/2) & (sz <= s.z_size);
    sx = sx(inside_idx);       
    sy = sy(inside_idx);
    sz = sz(inside_idx);
    amp = amp(inside_idx);

    % 若 scatterers 不足則補進新 speckle（補在 sy = -zipper_height/2 附近）
    num_new = s.num_scatterers - numel(sx);
    if num_new > 0
        new_sx = (rand(num_new,1) - 0.5) * s.x_size;
        new_sy = -s.zipper_height/2 + 0.1 * rand(num_new,1);  % 稍微進入畫面內的底部
        new_sz = rand(num_new,1) * s.z_size;

        sx = [sx; new_sx];
        sy = [sy; new_sy];
        sz = [sz; new_sz];
        new_amp = sqrt(-2 * log(rand(num_new,1)));
        amp = [amp; new_amp];
    end
end
function [original_image, moved_images] = generateAllImages(sx, sy, sz, amp, tx_list_left, rx_list_left, tx_list_right, rx_list_right, settings, X, Z, zipper_data)
    image_h = size(Z, 1);
    image_w = size(Z, 2);
    moved_images = cell(1, settings.num_iterations);

    % 左元件收 Frame 1 (t=0)
    fprintf('simulate origin image（左元件）\n');
    img_left = generateZipperImage(sx, sy, sz, amp, tx_list_left, rx_list_left, settings, X, Z, zipper_data);
    original_image = imresize(img_left, [image_h, image_w]);
    
    % 右元件收 Frame 2 (t=0)
    fprintf('simulate moved image 1（右元件）\n');
    img_right_t0 = generateZipperImage(sx, sy, sz, amp, tx_list_right, rx_list_right, settings, X, Z, zipper_data);
    moved_images{1} = imresize(img_right_t0, [image_h, image_w]);
    % 右元件收後續幀（模擬 speckle 往 y 軸移動）
    for iter = 2:settings.num_iterations
        fprintf('simulate moved image %d（右元件）\n', iter);
        [sx, sy, sz, amp] = updateScatterers(sx, sy, sz, amp, settings, iter - 1);
        img_right = generateZipperImage(sx, sy, sz, amp, tx_list_right, rx_list_right, settings, X, Z, zipper_data);
        moved_images{iter} = imresize(img_right, [image_h, image_w]);
    end
end
function img = generateZipperImage(sx, sy, sz, amp, tx_list, rx_list, s, X, Z, zipper_data)
    sigma_x = 0.3; sigma_z = 0.3; sigma_y0 = 0.2; alpha = 0.05;
    z_size = max(Z(:)); z_focal = z_size / 2;
    sigma_y = sigma_y0 + alpha*(Z - z_focal).^4;
    % 從 data 擷取發射/接收元件中心點
    num_elements = size(zipper_data, 1);  % 一共有幾個元件
    ele_center = zeros(num_elements, 3);
    for i = 1:num_elements
        ele_center(i, :) = zipper_data(i, 2:4);  % (x1, y1, z1)
    end

    flowField = zeros(size(X));
    for tx = tx_list
        active_element_y_center = ele_center(tx, 3);
        for rx = rx_list
            for i = 1:length(sx)
                tx_delay = norm([sx(i) sy(i) sz(i)] - ele_center(tx,:)) / s.c;
                rx_delay = norm([sx(i) sy(i) sz(i)] - ele_center(rx,:)) / s.c;
                total_delay = tx_delay + rx_delay;
                dir_tx = (sz(i) - ele_center(tx,3)) / norm([sx(i) - ele_center(tx,1), sy(i) - ele_center(tx,2), sz(i) - ele_center(tx,3)]);
                dir_rx = (sz(i) - ele_center(rx,3)) / norm([sx(i) - ele_center(rx,1), sy(i) - ele_center(rx,2), sz(i) - ele_center(rx,3)]);
                directivity = (dir_tx * dir_rx)^1.0;
                depth_idx = min(max(round((sz(i) / z_size) * size(sigma_y,1)),1), size(sigma_y,1));
                % elev_weight = exp(-(sy(i)^2) ./ (2 * sigma_y(depth_idx,1)^2));
                elev_weight = exp(-((sy(i) - active_element_y_center)^2) ./ (2 * sigma_y(depth_idx,1)^2));
                
                phase = 2*pi*s.fx*total_delay;
                psf = amp(i) * directivity * elev_weight .* exp(-(((X - sx(i)).^2)/(2*sigma_x^2) + ((Z - sz(i)).^2)/(2*sigma_z^2))) .* cos(2*pi*s.fx/s.c * (Z - sz(i)) + phase);
                % psf = amp(i) * directivity * elev_weight .* exp(-(((X - sx(i)).^2)/(2*sigma_x^2) + ((Z - sz(i)).^2)/(2*sigma_z^2)));
                flowField = flowField + psf;
            end
        end
    end

    envelope = abs(hilbert(flowField));
    envelope_norm = envelope / max(envelope(:));
    envelope_db = 20 * log10(envelope_norm); % 對數轉換＋閾值處理 
    img = max(envelope_db, -50);
end
function [disp_vec, cc_max, matched_block_pos] = calculateMatchedBlock_3D(orig_img, moved_img, k_block, window_w, window_h, s, iter, x, z)
    x0 = k_block(1); z0 = k_block(2);
    template = orig_img(z0:z0+window_h-1, x0:x0+window_w-1);
    
    if std(template(:)) == 0
        cc_max = 0; disp_vec = [0;0;0]; matched_block_pos = [x(x0), z(z0)]; return;
    end

    dx = s.dx; dz = s.dz;
    [h_img, w_img] = size(moved_img);

    pad_x = round((((window_w * 1.5)-3) / dx) / 2);
    pad_z = round((((window_h * 6)-3) / dz) / 2);
    pad_x = min([pad_x, x0 - 1, w_img - (x0 + window_w - 1)]);
    pad_z = min([pad_z, z0 - 1, h_img - (z0 + window_h - 1)]);
    if pad_x < 1 || pad_z < 1
        cc_max = 0; disp_vec = [0;0;0]; matched_block_pos = [x(x0), z(z0)]; return;
    end

    xmin = x0 - pad_x; xmax = x0 + window_w - 1 + pad_x;
    zmin = z0 - pad_z; zmax = z0 + window_h - 1 + pad_z;
    search_region = moved_img(zmin:zmax, xmin:xmax);

    C = normxcorr2(template, search_region);
    [cc_max, imax] = max(abs(C(:)));
    [ypeak, xpeak] = ind2sub(size(C), imax);
    yoffSet = ypeak - window_h;
    xoffSet = xpeak - window_w;

    matched_x0 = xmin + xoffSet;
    matched_z0 = zmin + yoffSet;

    matched_x0 = max(1, min(matched_x0, size(moved_img, 2) - window_w + 1));
    matched_z0 = max(1, min(matched_z0, size(moved_img, 1) - window_h + 1));

    dx_mm = (matched_x0 - x0) * dx;
    dz_mm = (matched_z0 - z0) * dz;
    dy_mm = (iter - 1) * s.Vmax * s.dt;

    disp_vec = [0; dy_mm; 0];
    matched_block_pos = [x(matched_x0), z(matched_z0)];
    
    % debug
    % 1) pad 與 search size
    fprintf('Frame %d: pad_x=%d, pad_z=%d, search=%dx%d px\n', ...
            iter, pad_x, pad_z, zmax-zmin+1, xmax-xmin+1);
    
    % 2) offset
    fprintf('  peak offset dx=%d dz=%d px\n', xoffSet, yoffSet);
    
    % 3) region 統計
    region = moved_img( matched_z0 : matched_z0 + window_h - 1, ...
                       matched_x0 : matched_x0 + window_w - 1 );
    fprintf('  std(T)=%.3f std(R)=%.3f mean(T)=%.3f mean(R)=%.3f\n', ...
            std(template(:)), std(region(:)), mean(template(:)), mean(region(:)));
    
    % 4) SNR
    Cabs = abs(C);
    mu = mean(Cabs(:)); sigma = std(Cabs(:));
    fprintf('  corr peak=%.3f, avg=%.3f, std=%.3f, SNR=%.2f\n', cc_max, mu, sigma, (cc_max-mu)/sigma);

end