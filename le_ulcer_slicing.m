%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function le_ulcer_slicing()

clc;

% parameters
bCutImage = 1;
nImageHalfWidth = 350; % in pixels | 350 (default) | stable parameter
aFilter = 'median'; % 'median' (default), 'lowpass', 'erode'

dGB = 10; % 10 (default), green - blue difference | detect ulcer
d2RGB = 80; % 80 (default), 2 * red - green - blue | exclude dark background

% get path
aPath = support_get_path();
aSubpath = support_fname({aPath, 'leprosy', 'TABLE_Aranz_Image'});

% get subjects
tSubjects = support_get_subjects(aSubpath);
nSubjects = length(tSubjects);

% smoothing filter
[fb, fa] = butter(4, 0.01, 'low');

% LP filter
[lfb, lfa] = butter(4, 0.02, 'low');

% loop subjects
for iSubject = 1:nSubjects
  aSubject = tSubjects{iSubject};

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

  % open figure
  hFigure = figure; 
  set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 2.0, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 

  % loop files
  nFiles = length(tFiles);
  pUlcerSize = zeros(1, nFiles);
  pDays = zeros(1, nFiles);
  nDateRef = datetime('now');
  for iFile = 1:nFiles
    aFile = tFiles{iFile};
    aFilename = support_fname({aSubpath, aSubject, aFile});

    % status
    fprintf(1, '%s: %s\n', aSubject, aFile);

    % get date
    aDate = aFile((end - 11):(end - 4));
    nDate = datetime([aDate(1:4), '-', aDate(5:6), '-', aDate(7:8)]);
    if iFile == 1
      nDateRef = nDate; 
    end
    nDateDif = days(nDate - nDateRef);
    pDays(iFile) = nDateDif;

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
    if strcmp(aFilter, 'median')
      D = [32, 32];
      J_R = medfilt2(I(:, :, 1), D);
      J_G = medfilt2(I(:, :, 2), D);
      J_B = medfilt2(I(:, :, 3), D);
      J = cat(3, J_R, J_G, J_B);
      H = double(J);
      S = double(J);
    elseif strcmp(aFilter, 'lowpass')
      J_R = I(:, :, 1);
      J_G = I(:, :, 2);
      J_B = I(:, :, 3);
      % smooth rows
      J_R = uint8(filtfilt(lfb, lfa, double(J_R)')');
      J_G = uint8(filtfilt(lfb, lfa, double(J_G)')');
      J_B = uint8(filtfilt(lfb, lfa, double(J_B)')');
      % smooth cols
      J_R = uint8(filtfilt(lfb, lfa, double(J_R)));
      J_G = uint8(filtfilt(lfb, lfa, double(J_G)));
      J_B = uint8(filtfilt(lfb, lfa, double(J_B)));
      % init
      J = cat(3, J_R, J_G, J_B);
      H = double(J);
      S = double(J);
    elseif strcmp(aFilter, 'erode')
      SE = strel('disk', 10); 
      J_R = imerode(I(:, :, 1), SE);
      J_G = imerode(I(:, :, 2), SE);
      J_B = imerode(I(:, :, 3), SE);
      J = cat(3, J_R, J_G, J_B);
      H = double(J);
      S = double(J);
    end

    H = abs(H(:, :, 2) - H(:, :, 3)) < dGB & 2 * H(:, :, 1) - H(:, :, 2) - H(:, :, 3) > d2RGB;

    % circle limit
    bCircleLimit = 1; % must be 1 always
    if bCircleLimit == 1
      x = sum(H, 1); 
      y = sum(H, 2); 
      x = filtfilt(fb, fa, x);
      y = filtfilt(fb, fa, y);
      cx = 0;
      cy = 0;
      bDebug = 0;
      if bDebug == 1
        figure;
        imshow(H); hold on;
        plot(1:size(H, 2), x + nImageHalfWidth, 'y');
        plot(y + nImageHalfWidth, 1:size(H, 1), 'c');
        plot(1:size(H, 2), iy * ones(1, size(H, 2)), 'y');
        plot(ix * ones(1, size(H, 1)), 1:size(H, 1), 'c');
        % borders
        % plot(1:size(H, 2), ones(size(x)) * (nImageHalfWidth - 50), 'y');
        % plot(1:size(H, 2), ones(size(x)) * (nImageHalfWidth + 50), 'y');
        % plot(ones(size(x)) * (nImageHalfWidth - 50), 1:size(H, 1), 'c');
        % plot(ones(size(x)) * (nImageHalfWidth + 50), 1:size(H, 1), 'c');
      end
      % range circles
      pR = 5:5:nImageHalfWidth;
      nR = length(pR);
      S = zeros(nR, 1);
      for iR = 1:nR
        R = pR(iR);
        s = sqrt(((-nImageHalfWidth:nImageHalfWidth) - cx) .^ 2 + ((-nImageHalfWidth:nImageHalfWidth)' - cy) .^ 2) < R;
        s = H .* s;
        S(iR) = sum(s(:));
      end
      xR = 10; % threshold
      dS = [0; diff(S)];
      hS = dS > xR;
      i = find(hS > 0, 1, 'first');
      if i < 20 % threshold
        iR = find(hS(i:end) < 1, 1, 'first') + i - 1;
        if isempty(iR)
          iR = length(pR);
        end
      else
        iR = 1;
      end
      R = pR(iR);
  
      bAddCircleToImage = 0; % 0 (default), 1 (visualisation only)
      if bAddCircleToImage == 1
        M = sqrt(((-nImageHalfWidth:nImageHalfWidth) - cx) .^ 2 + ((-nImageHalfWidth:nImageHalfWidth)' - cy) .^ 2) < R;
        C = sqrt(((-nImageHalfWidth:nImageHalfWidth) - cx) .^ 2 + ((-nImageHalfWidth:nImageHalfWidth)' - cy) .^ 2) > (R - 5);
        H = M .* H + C & M;
      else
        M = sqrt(((-nImageHalfWidth:nImageHalfWidth) - cx) .^ 2 + ((-nImageHalfWidth:nImageHalfWidth)' - cy) .^ 2) < R;
        H = M .* H;
      end
    end

    % dilate image
    bDilate = 0;
    if bDilate == 1
      nDilateDim = 8;
      SE = strel('disk', nDilateDim); 
      H = imdilate(H, SE); 
    end

    pUlcerSize(iFile) = sum(H(:)) / (nImageHalfWidth .^ 2);

    if iFile < 20
      subplot(5, 8, (iFile - 1) * 2 + 1); imshow(I);
      title(sprintf('day: %d', pDays(iFile)), 'FontWeight', 'normal', 'FontSize', 8);
      subplot(5, 8, (iFile - 1) * 2 + 2); imshow(H);
      title(sprintf('size: %1.4f', pUlcerSize(iFile)), 'FontWeight', 'normal', 'FontSize', 8);
    end
  end
  subplot(5, 8, [39, 40]); 
  x = pDays(pDays < 60);
  y = pUlcerSize(pDays < 60);
  p = fit(x(:), y(:), 'exp1');
  plot(p, x, y, '*'); box off;
  xlabel('size'); ylabel('days');
  set(gca, 'FontSize', 8);
  legend('off');
  title(sprintf('b = %1.3f', p.b), 'FontWeight', 'normal', 'FontSize', 8);
  
  % save image
  aDir = support_fname({aPath, 'leprosy', '_analysis', ['slicing_', aFilter]});
  if ~exist(aDir, 'dir')
    mkdir(aDir);
  end
  aFilename = support_fname({aDir, [aSubject, '.png']});
  print(hFigure, aFilename, '-dpng', '-r300');
  close(hFigure);  
end

end % end

%-------------------------------------------------------------------------------




