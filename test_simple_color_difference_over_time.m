%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_simple_color_difference_over_time()

clc;

% parameters
bCutImage = 1;
nImageHalfWidth = 350; % in pixels

% get path
aPath = support_get_path();
aSubpath = support_fname({aPath, 'leprosy', 'TABLE_Aranz_Image'});

% get subjects
tSubjects = support_get_subjects(aSubpath);
nSubjects = length(tSubjects);

% smoothing filter
% [fb, fa] = butter(2, 0.05, 'low');

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

  % create mask
  nFiles = length(tFiles);
  MASK = zeros(2 * nImageHalfWidth + 1, 2 * nImageHalfWidth + 1, nFiles);
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
    % median filter
    GF = 1.0 * I(:, :, 1) - I(:, :, 2) - I(:, :, 3); % fit this model
    GF = medfilt2(GF, [8, 8]);
    % threshold
    h = 32; % CHECK this
    Q = uint8((GF > h) * 255);
    % init
    MASK(:, :, iFile) = double(GF > h);
  end
  MASK = sum(MASK, 3);
  MASK = MASK > 1.0; % MASK threshold (?)
  % fit ellipse
  MASK = fit_ellipse(MASK);

  % static MASK
  bStaticMASK = 0;
  if bStaticMASK == 1
    R = 200;
    MASK = uint8(sqrt((-nImageHalfWidth:nImageHalfWidth) .^ 2 + (-nImageHalfWidth:nImageHalfWidth)' .^ 2) < R);
  end

  % loop files
  nFiles = length(tFiles);
  pUlcerSize = zeros(nFiles, 1);
  for iFile = 1:12 % 1:nFiles
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

    % color difference
    % G = I(:, :, 1) - I(:, :, 2) - I(:, :, 3); 
    % G = repmat(G, 1, 1, 3);

    GF = 1.2 * I(:, :, 1) - I(:, :, 2) - I(:, :, 3); % fit this model
    GF = medfilt2(GF, [8, 8]);
 
    % threshold
    h = 32;
    Q = uint8((GF > h) * 255);

    Q_MASK = Q .* uint8(MASK);

    % ulcer size
    nUlcerSize = sum(Q_MASK(:) > 0) / length(Q_MASK(:));
    pUlcerSize(iFile) = nUlcerSize;

    subplot(4, 4, iFile);
    imshow(Q_MASK); % Q, I

    hold on;

    plot(1:(2 * nImageHalfWidth), zeros(2 * nImageHalfWidth, 1) + nImageHalfWidth, 'LineWidth', 1, 'Color', 'c');
    plot(zeros(2 * nImageHalfWidth, 1) + nImageHalfWidth, 1:(2 * nImageHalfWidth), 'LineWidth', 1, 'Color', 'c');

    title(sprintf('Ulcer: %1.4f', nUlcerSize), 'FontWeight', 'normal');

    % QI = uint8(Q > 0) .* I;
    % subplot(2, 4, 4);
    % imshow(QI);

  end
  subplot(4, 4, 13);
  imagesc(MASK);

  subplot(4, 4, [15, 16]);
  plot(pUlcerSize, '-*'); box off;
 
  return
end


end % end

%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function CI = fit_ellipse(I)

threshold = 0.08;

% init
[nHeight, nWidth] = size(I);
BI = I;
x = I;

% fitting ellipse 
pys = sum(x, 1); % x is image of ulcer
pxs = sum(x, 2);
pyg = double(pys > threshold * max(pys)) * size(x, 1);
pxg = double(pxs > threshold * max(pxs)) * size(x, 2);

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
bPlotMASK = 0;
if bPlotMASK == 1
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
bPlotMASK = 0;
if bPlotMASK == 1
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
bPlotMASK = 0;
if bPlotMASK == 1
  subplot(2, 2, 3); 
  imshow(CI); hold on;
  plot(pyg, 'Color', 'y', 'LineWidth', 0.5);
  plot(pxg, 1:length(pxg), 'Color', 'g', 'LineWidth', 0.5);
  title(sprintf('ovl = %1.2f%%', nOvlMax), 'FontWeight', 'normal');
end

end % end

%-------------------------------------------------------------------------------




