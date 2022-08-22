%-------------------------------------------------------------------------------
% Function
% LINK:
%   * https://uk.mathworks.com/discovery/image-segmentation.html
%   * https://uk.mathworks.com/help/vision/ug/segment-3d-brain-tumor-using-deep-learning.html
%-------------------------------------------------------------------------------
function le_test_load_images_eval()

clc;

% segmentatyion parameters
threshold = 0.08;

% get path
aPath = support_get_path();
aSubpath = support_fname({aPath, 'leprosy', 'TABLE_Aranz_Image'});

% get subjects
tSubjects = support_get_subjects(aSubpath);
nSubjects = length(tSubjects);

% loop subjects
for iSubject = 1:nSubjects
  aSubject = tSubjects{iSubject};
  tFiles = [];
  a = dir(support_fname({aSubpath, aSubject}));
  k = 1;
  for i = 1:length(a)
    if contains(a(i).name, '.png') || contains(a(i).name, '.jpg')
      tFiles{k} = a(i).name;
      k = k + 1;
    end
  end
  % loop files
  nFiles = length(tFiles);
  for iFile = 1:nFiles
    aFile = tFiles{iFile};
    aFilename = support_fname({aSubpath, aSubject, aFile});
    % load image
    I = imread(aFilename);
    nWidth = size(I, 2);
    nHeight = size(I, 1);

    % plot
    subplot(2, 4, 1); 
    imshow(I);
    title('Original image', 'FontWeight', 'normal');

    % convert RGB to CIE 1976 L*a*b
    he = I;
    lab_he = rgb2lab(he);

    % subplot(2, 4, 2); 
    % imshow(lab_he);

    numColors = 5;

    ab = lab_he(:, :, 2:3);
    ab = im2single(ab);

        % K-means clustering based image segmentation
    pixel_labels = imsegkmeans(ab, numColors, 'NumAttempts', 3);

    B2 = labeloverlay(he, pixel_labels);

    % plot
    subplot(2, 4, 2); 
    imshow(B2);
    title('Labeled Image a*b*', 'FontWeight', 'normal');

    % select cluster
    s_cols = zeros(1, numColors);
    s_rows = zeros(1, numColors);
    for k = 1:numColors
      mask = pixel_labels == k;
      % quantify clusters
      cols = sum(mask, 1);
      rows = sum(mask, 2)';
      % threshold
      th_cols = 0.05 * max(cols) * ones(1, nWidth);
      th_rows = 0.05 * max(rows) * ones(1, nHeight);
      % above threshold
      s_cols(k) = mean(cols > th_cols);
      s_rows(k) = mean(rows > th_rows);
      % plot
      bPlot = 0; 
      if bPlot == 1 
        subplot(4, 4, 3 - 1 + k); 
        plot(cols); hold on; plot(rows); 
        plot(th_cols, 'r'); hold on; plot(th_rows, 'k');
      end
    end
    s = s_cols + s_rows;
    [~, j] = sort(s);

    % init
    ulcer = pixel_labels == j(1);
    feet = pixel_labels == j(1) | pixel_labels == j(2) | pixel_labels == j(3);

    % plot
    subplot(2, 4, 3); 
    imshow(ulcer);
    title('ulcer', 'FontWeight', 'normal');

    subplot(2, 4, 4); 
    imshow(feet);
    title('feet', 'FontWeight', 'normal');

    % median filtering
    x = ulcer;
    x = medfilt2(x, [32, 32]);

    % plot
    subplot(2, 4, 5); 
    imshow(x);
    title('Median filter images', 'FontWeight', 'normal');

    % deliniate cluster
    x = double(x > 0);
    BI = x;

    % plot
    subplot(2, 4, 6); 
    imshow(x); hold on;
    py = sum(x, 1) + size(x, 1) / 2;
    px = sum(x, 2) + size(x, 2) / 2;
    plot(py, 'Color', 'y', 'LineWidth', 0.5);
    plot(px, 1:size(x, 1), 'Color', 'g', 'LineWidth', 0.5);
    title('Cluster borders', 'FontWeight', 'normal');

    % fitting ellipse 
    subplot(2, 4, 7); 
    imshow(x); hold on;

    pys = sum(x, 1);
    pxs = sum(x, 2);
    pyg = double(pys > threshold * max(pys)) * size(x, 1);
    pxg = double(pxs > threshold * max(pxs)) * size(x, 2);

    pyg_min = find(pyg > 0, 1, 'first');
    pyg_max = find(pyg > 0, 1, 'last');

    pxg_min = find(pxg > 0, 1, 'first');
    pxg_max = find(pxg > 0, 1, 'last');

    plot(pyg, 'Color', 'y', 'LineWidth', 0.5);
    plot(pxg, 1:size(x, 1), 'Color', 'g', 'LineWidth', 0.5);

    % [~, iy] = max(py);
    % [~, ix] = max(px);

    % evaluate eccentricity 
    ECC = 99;
    pOvls = zeros(1, ECC);
    pEccs = zeros(1, ECC);
    for ecc = 1:ECC

      eccentricity = ecc * 0.01;
      pEccs(ecc) = eccentricity;

      % fitting elipse
      x1 = pyg_min;
      x2 = pyg_max;
      y1 = pxg_min;
      y2 = pxg_max;
  
      a = (1/2) * sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
      b = a * sqrt(1 - eccentricity ^ 2);
      t = linspace(0, 2 * pi, 1000); 
      % compute angles relative to (x1, y1)
      angles = atan2(y2 - y1, x2 - x1);
      x = (x1 + x2) / 2 + a * cos(t) * cos(angles) - b * sin(t) * sin(angles);
      y = (y1 + y2) / 2 + a * cos(t) * sin(angles) + b * sin(t) * cos(angles);

      % hide now
      %%%% plot(x, y, 'Color', 'r', 'LineWidth', 2);

      % make image circle
      CI = zeros(nHeight, nWidth);
  
      gx = round(x);
      gy = round(y);
  
      for i = 1:length(t)
        CI(gy(i), ((gx(i) - 1):(gx(i) + 1))) = 1;
      end
  
      gy_min = min(gy);
      gy_max = max(gy);
      for i = gy_min:gy_max
        K = find(CI(i, :) > 0);
        if ~isempty(K)
          CI(i, K(1):K(end)) = 1;
        end
      end

      % hide now
      % subplot(2, 4, 7); 
      % imshow(CI); hold on;
      % plot(pyg, 'Color', 'y', 'LineWidth', 0.5);
      % plot(pxg, 1:length(pxg), 'Color', 'g', 'LineWidth', 0.5);

      % overlap
      g_and = BI & CI;
      g_or = BI | CI;
      ovl = sum(g_and(:)) / sum(g_or(:));
      pOvls(ecc) = ovl;
    end

    % get optimal eccentricity
    [~, j] = max(pOvls);
    eccentricity = pEccs(j);

    % plot
    subplot(2, 4, 8); 
    plot(pEccs, pOvls, 'k.'); box off;
    title(sprintf('eccentricity = %1.2f', eccentricity), 'FontWeight', 'normal');

    % fitting elipse
    x1 = pyg_min;
    x2 = pyg_max;
    y1 = pxg_min;
    y2 = pxg_max;

    a = (1/2) * sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
    b = a * sqrt(1 - eccentricity ^ 2);
    t = linspace(0, 2 * pi, 1000); 
    % compute angles relative to (x1, y1)
    angles = atan2(y2 - y1, x2 - x1);
    x = (x1 + x2) / 2 + a * cos(t) * cos(angles) - b * sin(t) * sin(angles);
    y = (y1 + y2) / 2 + a * cos(t) * sin(angles) + b * sin(t) * cos(angles);

    % plot
    subplot(2, 4, 6); 
    plot(x, y, 'Color', 'r', 'LineWidth', 2);

    % make image circle
    CI = zeros(nHeight, nWidth);

    gx = round(x);
    gy = round(y);

    for i = 1:length(t)
      CI(gy(i), ((gx(i) - 1):(gx(i) + 1))) = 1;
    end

    gy_min = min(gy);
    gy_max = max(gy);
    for i = gy_min:gy_max
      K = find(CI(i, :) > 0);
      if ~isempty(K)
        CI(i, K(1):K(end)) = 1;
      end
    end

    % plot
    subplot(2, 4, 7); 
    imshow(CI); hold on;
    plot(pyg, 'Color', 'y', 'LineWidth', 0.5);
    plot(pxg, 1:length(pxg), 'Color', 'g', 'LineWidth', 0.5);

    % plot | overlay ellipse and original image
    subplot(2, 4, 1); hold on;
    plot(x, y, 'Color', 'r', 'LineWidth', 3);

    return

  end

end


end % end

%-------------------------------------------------------------------------------




