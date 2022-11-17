%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_simple_color_difference_model_fitting()

clc;

% parameters
bCutImage = 1;
nImageHalfWidth = 350; % in pixels | 350 (default) | stable parameter

nRedWeight = 1.15; % 1.15 (default) | stable parameter

MASK_nRedWeight = 1.0; % 1.0 (default) | stable parameter
MASK_nThreshold = 0.0; % 1.0 (default) | sensitive parameter | if decreasing this parameter then consider increasing BLOB_nThreshold

BLOB_nThreshold = 0.005; % 0.005 (default) | roughly, ratio of white pixels between blobs norm to max

MASK_pEllipseRatio = [2, 3]; % [2, 3] (default) | stable parameter
MEDF_nSize = [8, 8]; % [8, 8] (default) | stable parameter
MEDF_nThreshold = 32; % 32 (default) | stable parameter

% get path
aPath = support_get_path();
aSubpath = support_fname({aPath, 'leprosy', 'TABLE_Aranz_Image'});

% get subjects
tSubjects = support_get_subjects(aSubpath);
nSubjects = length(tSubjects);

% smoothing filter
[fb, fa] = butter(4, 0.05, 'low');

% loop subjects
for iSubject = 1:1 % nSubjects
  aSubject = tSubjects{iSubject};
  % debug
  % aSubject = 'S-283';

  % status
  fprintf(1, '%s\n', aSubject);

  % open figure
  % hFigure = figure; 
  % set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 2.0, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 

  % get files
  tFiles = [];
  a = dir(support_fname({aSubpath, aSubject}));
  k = 1;
  for i = 1:length(a)
    if contains(a(i).name, '.png') || contains(a(i).name, '.jpg')
      tFiles{k} = a(i).name;
      k = k + 1;
    end
  end

  % get initial image
  aFile = tFiles{1};
  aFilename = support_fname({aSubpath, aSubject, aFile});
  % load image
  I = imread(aFilename);
  nWidth = size(I, 2);
  nHeight = size(I, 1);
	% cut image
  if bCutImage == 1
    d = nImageHalfWidth;
		x = nWidth / 2;
		y = nHeight / 2;
		I = I(:, :, :); 
		I = I((y - d):(y + d), (x - d):(x + d), :);
  end


  % model

  GG = 0; % previous
 

  gW = 0.5:0.025:1.1;
  gR = zeros(1, length(gW));
  gG = zeros(1, length(gW));

  figure(1);

  for iW = 1:length(gW) 
    subplot(4, 7, iW);
    J = gW(iW) * I(:, :, 1) - I(:, :, 2) - I(:, :, 3);
    J = double(J > 0);

    x = sum(J, 1); 
    y = sum(J, 2); 

    x = filtfilt(fb, fa, x);
    y = filtfilt(fb, fa, y);

    % min to max
    h = 5 / 2;

    ix0 = find(x > h, 1, 'first');
    ix1 = find(x > h, 1, 'last');

    iy0 = find(y > h, 1, 'first');
    iy1 = find(y > h, 1, 'last');

    ix = (ix1 - ix0) / 2 + ix0;
    iy = (iy1 - iy0) / 2 + iy0;

    % max
%     [~, ix] = max(x);
%     [~, iy] = max(y);

    cx = ix - nImageHalfWidth;
    cy = iy - nImageHalfWidth;
    
    pR = 5:5:nImageHalfWidth;
    nR = length(pR);
    S = zeros(nR, 1);
    for iR = 1:nR
      R = pR(iR);
      s = sqrt(((-nImageHalfWidth:nImageHalfWidth) - cx) .^ 2 + ((-nImageHalfWidth:nImageHalfWidth)' - cy) .^ 2) < R;
      s = J .* s;
      S(iR) = sum(s(:));
      % subplot(4, 7, i);
      % imshow(s);
    end
    xR = 5;
    dS = [0; diff(S)];
    [~, i] = max(dS);

    iR = find(dS(i:end) < xR, 1, 'first') + i;

    R = pR(iR);

    s = sqrt(((-nImageHalfWidth:nImageHalfWidth) - cx) .^ 2 + ((-nImageHalfWidth:nImageHalfWidth)' - cy) .^ 2) < R;

    % bPlotA1 = 0;
    % if bPlotA1 == 1
    %   subplot(2, 2, 1); imshow(s); hold on;
    %   plot(1:size(J, 2), iy * ones(1, size(J, 2)), 'y');
    %   plot(ix * ones(1, size(J, 1)), 1:size(J, 1), 'c');
    % 
    %   subplot(2, 2, 2); imshow(J); hold on;
    %   plot(1:size(J, 2), iy * ones(1, size(J, 2)), 'y');
    %   plot(ix * ones(1, size(J, 1)), 1:size(J, 1), 'c');
    % 
    %   subplot(2, 2, 3); imshow(s - J); hold on;
    %   plot(1:size(J, 2), iy * ones(1, size(J, 2)), 'y');
    %   plot(ix * ones(1, size(J, 1)), 1:size(J, 1), 'c');
    % end

    imshow(s - J);

    G = sum(J(:)) / R .^ 2;

    gG(iW) = G;
    gR(iW) = R;

    bBREAK = 0;
    if (G - GG) < 0 && G / GG < 0.95
      bBREAK = 1;
    end

    GG = G;
    
    title(sprintf('%1.2f | %d | %1.2f', gW(iW), R, G), 'FontWeight', 'normal');

    if bBREAK == 1
      iW = iW + 1;
      break
    end

  end
  subplot(4, 7, iW);
  plot(gW, gG, '-*');
  title('G', 'FontWeight', 'normal');
  subplot(4, 7, iW + 1);
  plot(gW, gR, '-*');
  title('R', 'FontWeight', 'normal');

  % We have here: W(red) and R(image_1) for MASK (!)
  % Then we apply W(red) for second image, detect center of image and apply MASK with R(image_1)

  % init R
  R = gR(iW - 2);
  W = gW(iW - 2);


  figure(2);

  % loop files
  nFiles = length(tFiles);
  pUlcerSize = zeros(nFiles, 1);
  for iFile = 1:nFiles
    aFile = tFiles{iFile};
    aFilename = support_fname({aSubpath, aSubject, aFile});
    % load image
    I = imread(aFilename);
    nWidth = size(I, 2);
    nHeight = size(I, 1);
		
		% cut image
    if bCutImage == 1
      d = nImageHalfWidth;
			x = nWidth / 2;
			y = nHeight / 2;
			I = I(:, :, :); 
			I = I((y - d):(y + d), (x - d):(x + d), :);
    end

    J = W * I(:, :, 1) - I(:, :, 2) - I(:, :, 3);
    J = double(J > 0);

    x = sum(J, 1); 
    y = sum(J, 2); 

    x = filtfilt(fb, fa, x);
    y = filtfilt(fb, fa, y);

    % min to max
    h = 5 / 2;

    ix0 = find(x > h, 1, 'first');
    ix1 = find(x > h, 1, 'last');

    iy0 = find(y > h, 1, 'first');
    iy1 = find(y > h, 1, 'last');

    ix = (ix1 - ix0) / 2 + ix0;
    iy = (iy1 - iy0) / 2 + iy0;

    % max
%     [~, ix] = max(x);
%     [~, iy] = max(y);

    cx = ix - nImageHalfWidth;
    cy = iy - nImageHalfWidth;

    s = sqrt(((-nImageHalfWidth:nImageHalfWidth) - cx) .^ 2 + ((-nImageHalfWidth:nImageHalfWidth)' - cy) .^ 2) < R;
    MASK = s;




    % ulcer size
% %     nUlcerSize = sum(Q_MASK(:) > 0) / length(Q_MASK(:));
% %     pUlcerSize(iFile) = nUlcerSize;

    subplot(4, 8, (iFile - 1) * 2 + 1);
    imshow(I); 

    subplot(4, 8, (iFile - 1) * 2 + 2);
    imshow(J .* MASK); % imshow(MASK); 


    hold on;
% % % % 
% % % %     % cross hair
% % % %     bCrossHair = 0;
% % % %     if bCrossHair == 1
% % % %       plot(1:(2 * nImageHalfWidth), zeros(2 * nImageHalfWidth, 1) + nImageHalfWidth, 'LineWidth', 1, 'Color', 'c');
% % % %       plot(zeros(2 * nImageHalfWidth, 1) + nImageHalfWidth, 1:(2 * nImageHalfWidth), 'LineWidth', 1, 'Color', 'c');
% % % %       title(sprintf('Ulcer: %1.4f', nUlcerSize), 'FontWeight', 'normal');
% % % %     end

    % break
    if iFile > 14
      break
    end
  end
  subplot(4, 8, 31);
  imshow(MASK);
  subplot(4, 8, 32);
  plot(pUlcerSize, '-*'); box off;
 
  % save figure
  aFilename = support_fname({aPath, 'leprosy', '_analysis', 'simple_color_difference', [aSubject, '.png']});
  print(hFigure, aFilename, '-dpng', '-r300');
  close(hFigure); 
end

end % end

%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function CI = fit_ellipse_eccentricity(I, bDebug)

% init
[nHeight, nWidth] = size(I);
BI = I;
x = I;

% fitting ellipse 
pys = sum(x, 1); % x is image of ulcer
pxs = sum(x, 2);

% version 1.0
% threshold = 0.08; 
% pyg = double(pys > threshold * max(pys)) * size(x, 1);
% pxg = double(pxs > threshold * max(pxs)) * size(x, 2);

% version 1.1 | threshold 1.0
pyg = zeros(size(pys)); [~, i] = max(pys); jR = find(pys(i:end) < 1.0, 1, 'first'); jL = find(pys(i:-1:1) < 1.0, 1, 'first'); pyg((i - jL + 1):(i + jR - 1)) = 1;
pxg = zeros(size(pxs)); [~, i] = max(pxs); jR = find(pxs(i:end) < 1.0, 1, 'first'); jL = find(pxs(i:-1:1) < 1.0, 1, 'first'); pxg((i - jL + 1):(i + jR - 1)) = 1;

pyg_min = find(pyg > 0, 1, 'first');
pyg_max = find(pyg > 0, 1, 'last');
pxg_min = find(pxg > 0, 1, 'first');
pxg_max = find(pxg > 0, 1, 'last');

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

  % overlap
  g_and = BI & CI;
  g_or = BI | CI;
  ovl = sum(g_and(:)) / sum(g_or(:));
  pOvls(ecc) = ovl;
end

% get optimal eccentricity
[nOvlMax, j] = max(pOvls);
eccentricity = pEccs(j);

% plot
if bDebug == 1
  figure; 
  subplot(2, 2, 2); 
  plot(pEccs, pOvls, 'k.'); hold on;
  plot(pEccs(j), pOvls(j), 'LineWidth', 2, 'Color', 'r', 'Marker', 'o'); box off;
  xlabel('eccentricity'); ylabel('overlap');
  title(sprintf('ecc = %1.2f, ovl = %1.2f', eccentricity, nOvlMax), 'FontWeight', 'normal');
end

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

% plot | overlay ellipse and original image
if bDebug == 1
  subplot(2, 2, 1); 
  imshow(I); hold on;
  plot(x, y, 'Color', 'r', 'LineWidth', 3);
  title('MASK', 'FontWeight', 'normal');
end

% fitted mask
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
if bDebug == 1
  subplot(2, 2, 3); 
  imshow(CI); hold on;
  plot(pyg, 'Color', 'y', 'LineWidth', 0.5);
  plot(pxg, 1:length(pxg), 'Color', 'g', 'LineWidth', 0.5);
  title(sprintf('ovl = %1.2f%%', nOvlMax), 'FontWeight', 'normal');
end

end % end

%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function CI = fit_ellipse_angle(I, xy, bDebug)

% init
[nHeight, nWidth] = size(I);
BI = I;
x = I;

% fitting ellipse 
pys = sum(x, 1); % x is image of ulcer
pxs = sum(x, 2);

% version 1.0
% threshold = 0.08; 
% pyg = double(pys > threshold * max(pys)) * size(x, 1);
% pxg = double(pxs > threshold * max(pxs)) * size(x, 2);

% version 1.1 | threshold 1.0
pyg = zeros(size(pys)); [~, i] = max(pys); jR = find(pys(i:end) < 1.0, 1, 'first'); jL = find(pys(i:-1:1) < 1.0, 1, 'first'); pyg((i - jL + 1):(i + jR - 1)) = 1;
pxg = zeros(size(pxs)); [~, i] = max(pxs); jR = find(pxs(i:end) < 1.0, 1, 'first'); jL = find(pxs(i:-1:1) < 1.0, 1, 'first'); pxg((i - jL + 1):(i + jR - 1)) = 1;

pyg_min = find(pyg > 0, 1, 'first');
pyg_max = find(pyg > 0, 1, 'last');
pxg_min = find(pxg > 0, 1, 'first');
pxg_max = find(pxg > 0, 1, 'last');

% evaluate angle 
ECC = 90; % 180
pOvls = zeros(1, ECC);
pPhis = zeros(1, ECC);
for ecc = 1:ECC

  phi = ecc ;
  pPhis(ecc) = phi;

  % fitting elipse
  x1 = pyg_min;
  x2 = pyg_max;
  y1 = pxg_min;
  y2 = pxg_max;

  rx = (x2 - x1) / xy(1);
  ry = (y2 - y1) / xy(2);

  t = linspace(0, 2 * pi, 1000); 
  phi_rad = (phi / 180) * pi;
  x = (x1 + x2) / 2 + rx * cos(t) * cos(phi_rad) - ry * sin(t) * sin(phi_rad);
  y = (y1 + y2) / 2 + rx * cos(t) * sin(phi_rad) + ry * sin(t) * cos(phi_rad);

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

  % overlap
  g_and = BI & CI;
  g_or = BI | CI;
  ovl = sum(g_and(:)) / sum(g_or(:));
  pOvls(ecc) = ovl;
end

% get optimal eccentricity
[nOvlMax, j] = max(pOvls);
phi = pPhis(j);

% plot
if bDebug == 1
  figure; 
  subplot(2, 2, 2); 
  plot(pPhis, pOvls, 'k.'); hold on;
  plot(pPhis(j), pOvls(j), 'LineWidth', 2, 'Color', 'r', 'Marker', 'o'); box off;
  xlabel('eccentricity'); ylabel('overlap');
  title(sprintf('angle = %1.2f, ovl = %1.2f', phi, nOvlMax), 'FontWeight', 'normal');
end

% fitting elipse
x1 = pyg_min;
x2 = pyg_max;
y1 = pxg_min;
y2 = pxg_max;

rx = (x2 - x1) / xy(1);
ry = (y2 - y1) / xy(2);

t = linspace(0, 2 * pi, 1000); 
phi_rad = (phi / 180) * pi;
x = (x1 + x2) / 2 + rx * cos(t) * cos(phi_rad) - ry * sin(t) * sin(phi_rad);
y = (y1 + y2) / 2 + rx * cos(t) * sin(phi_rad) + ry * sin(t) * cos(phi_rad);

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

% plot | overlay ellipse and original image
if bDebug == 1
  subplot(2, 2, 1); 
  imshow(I); hold on;
  plot(x, y, 'Color', 'r', 'LineWidth', 3);
  title('MASK', 'FontWeight', 'normal');
end

% fitted mask
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
if bDebug == 1
  subplot(2, 2, 3); 
  imshow(CI); hold on;
  plot(pyg, 'Color', 'y', 'LineWidth', 0.5);
  plot(pxg, 1:length(pxg), 'Color', 'g', 'LineWidth', 0.5);
  title(sprintf('ovl = %1.2f%%', nOvlMax), 'FontWeight', 'normal');
end

end % end

%-------------------------------------------------------------------------------




