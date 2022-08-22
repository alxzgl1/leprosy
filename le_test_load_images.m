%-------------------------------------------------------------------------------
% Function
% LINK:
%   * https://uk.mathworks.com/discovery/image-segmentation.html
%   * https://uk.mathworks.com/help/vision/ug/segment-3d-brain-tumor-using-deep-learning.html
%-------------------------------------------------------------------------------
function le_test_load_images()

clc;

threshold = 0.08;

eccentricity = 0.85;

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

    subplot(2, 4, 1); 
    imshow(I);

    he = I;
    lab_he = rgb2lab(he);

    % subplot(2, 4, 2); 
    % imshow(lab_he);

    numColors = 5;

    ab = lab_he(:, :, 2:3);
    ab = im2single(ab);
    pixel_labels = imsegkmeans(ab, numColors, "NumAttempts", 3);

    B2 = labeloverlay(he, pixel_labels);
    subplot(2, 4, 2); 
    imshow(B2);
    title("Labeled Image a*b*");


    mask4 = pixel_labels == 4;
    cluster1 = he .* uint8(mask4);
    subplot(2, 4, 3); 
    imshow(cluster1)
    title("Objects in Cluster 1");

%     L = lab_he(:, :, 1);
%     L_blue = L .* double(mask4);
%     L_blue = rescale(L_blue);
%     idx_light_blue = imbinarize(nonzeros(L_blue));
% 
%     blue_idx = find(mask4);
%     mask_dark_blue = mask4;
%     mask_dark_blue(blue_idx(idx_light_blue)) = 0;
%     
%     blue_nuclei = he .* uint8(mask_dark_blue);
%     subplot(2, 4, 5); 
%     imshow(blue_nuclei)


    x = rgb2gray(cluster1);
    x = medfilt2(x, [32, 32]);

    subplot(2, 4, 4); 
    imshow(x);

    x = double(x > 0);
    BI = x;

    subplot(2, 4, 5); 
    imshow(x); hold on;
    py = sum(x, 1) + size(x, 1) / 2;
    px = sum(x, 2) + size(x, 2) / 2;

    plot(py, 'Color', 'y', 'LineWidth', 0.5);
    plot(px, 1:size(x, 1), 'Color', 'g', 'LineWidth', 0.5);


    subplot(2, 4, 6); 
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

    

    % Make equations:
    x1 = pyg_min;
    x2 = pyg_max;
    y1 = pxg_min;
    y2 = pxg_max;

    a = (1/2) * sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
    b = a * sqrt(1 - eccentricity ^ 2);
    t = linspace(0, 2 * pi, 1000); 
    % Compute angles relative to (x1, y1).
    angles = atan2(y2 - y1, x2 - x1);
    x = (x1 + x2) / 2 + a * cos(t) * cos(angles) - b * sin(t) * sin(angles);
    y = (y1 + y2) / 2 + a * cos(t) * sin(angles) + b * sin(t) * cos(angles);

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
      CI(i, K(1):K(end)) = 1;
    end

    subplot(2, 4, 7); 
    imshow(CI); hold on;
    plot(pyg, 'Color', 'y', 'LineWidth', 0.5);
    plot(pxg, 1:length(pxg), 'Color', 'g', 'LineWidth', 0.5);

    % overlap
    g_and = BI & CI;
    g_or = BI | CI;
    ovl = sum(g_and(:)) / sum(g_or(:));

    title(sprintf('%1.4f', ovl), 'FontWeight', 'normal');


  subplot(2, 4, 1); hold on;
  plot(x, y, 'Color', 'r', 'LineWidth', 3);



    % gaussian fit
%     x = double(mean(blue_nuclei, 3));
%     GMModel = fitgmdist(x, 2);


    % grayscale
%     I = rgb2gray(I);

%     subplot(2, 2, 2); 
%     imshow(I);

    % I = imbinarize(I); 
    
    % active contour
    % mask = zeros(size(I));
    % mask(25:end-25, 25:end-25) = 1;
    % I = activecontour(I, mask);
    
    % multi threshold
%     level = multithresh(I);
%     seg_I = imquantize(I, level);
%     subplot(2, 2, 3); 
%     imshow(seg_I, []);

    % image
%     subplot(2, 2, 2); imshow(I(:, :, 1));
%     subplot(2, 2, 3); imshow(I(:, :, 2));
%     subplot(2, 2, 4); imshow(I(:, :, 3));

    return
    o = 0;
  end

end


end % end

%-------------------------------------------------------------------------------




