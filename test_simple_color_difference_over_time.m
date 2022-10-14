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
% [fb, fa] = butter(2, 0.2, 'low');

% loop subjects
for iSubject = 2:nSubjects
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

  % create mask by averaging subject's images 
  % parameters: (1) medfilt([8, 8]), (2) h = 32, (3) group MASK > 1
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
    % init
    MASK(:, :, iFile) = double(GF > 32); % arbitrary threshold of 32
  end
  MASK = sum(MASK, 3);
  MASK = MASK > 1.0; % threshold group MASK

  % exclude peripheral (artificial) blobs 
  % parameters: (1) diff(S) < 5 (0 = no white pixels between blobs) 
  bExcludePeripheralBlobs = 1;
  if bExcludePeripheralBlobs == 1
    bDebug = 0; 
    if bDebug == 1, figure; subplot(2, 2, 1); imshow(MASK); end
    pR = 5:5:nImageHalfWidth;
    nR = length(pR);
    S = zeros(nR, 1);
    for i = 1:nR
      R = pR(i);
      s = sqrt((-nImageHalfWidth:nImageHalfWidth) .^ 2 + (-nImageHalfWidth:nImageHalfWidth)' .^ 2) < R;
      s = MASK .* s;
      S(i) = sum(s(:));
    end
    % plot(diff(S))
    i = find(diff(S) < 5, 1, 'first'); % arbitrary threshold to separate central and peripheral (artificial) blobs 
    % make circle mask with optimal radius
    R = pR(i);
    s = sqrt((-nImageHalfWidth:nImageHalfWidth) .^ 2 + (-nImageHalfWidth:nImageHalfWidth)' .^ 2) < R;
    MASK = MASK .* s;
    if bDebug == 1, subplot(2, 2, 2); imshow(MASK); end
  end

  % fit ellipse
  bEllipseMASK = 1;
  if bEllipseMASK == 1
    bDebug = 0;
    MASK = fit_ellipse(MASK, bDebug);
  end

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

    GF = 1.15 * I(:, :, 1) - I(:, :, 2) - I(:, :, 3); % fit this model
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

    % cross hair
    plot(1:(2 * nImageHalfWidth), zeros(2 * nImageHalfWidth, 1) + nImageHalfWidth, 'LineWidth', 1, 'Color', 'c');
    plot(zeros(2 * nImageHalfWidth, 1) + nImageHalfWidth, 1:(2 * nImageHalfWidth), 'LineWidth', 1, 'Color', 'c');
    title(sprintf('Ulcer: %1.4f', nUlcerSize), 'FontWeight', 'normal');
  end
  subplot(4, 4, 13);
  imshow(MASK);
  subplot(4, 4, 14);
  plot(pUlcerSize, '-*'); box off;
 
  return
end


end % end

%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function CI = fit_ellipse(I, bDebug)

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




