%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_ulcer_slicing_plots()

clc;

% parameters
bCutImage = 1;
nImageHalfWidth = 350; % in pixels | 350 (default) | stable parameter

% get path
aPath = support_get_path();
aSubpath = support_fname({aPath, 'leprosy', 'TABLE_Aranz_Image'});

% get subjects
tSubjects = support_get_subjects(aSubpath);
nSubjects = length(tSubjects);

% smoothing filter
[fb, fa] = butter(4, 0.01, 'low');

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

  % loop files
  nDateRef = datetime('now');
  nFiles = length(tFiles);
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
    D = [20, 20];
    J_R = medfilt2(I(:, :, 1), D);
    J_G = medfilt2(I(:, :, 2), D);
    J_B = medfilt2(I(:, :, 3), D);
    J = cat(3, J_R, J_G, J_B);
    H = double(J);
    S = double(J);

    % open figure
    hFigure = figure; 
    set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 2.0, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 

    % decompose
    j = 1;
    for i = 100:50:((2 * nImageHalfWidth) - 50)
      subplot(4, 4, j + 4);
      bCols = 0;
      if bCols == 1
        RI = squeeze(S(:, i, 1));
        GI = squeeze(S(:, i, 2));
        BI = squeeze(S(:, i, 3));
        J(:, (i + 1):(i + 4), :) = 0;
      else
        RI = squeeze(S(i, :, 1));
        GI = squeeze(S(i, :, 2));
        BI = squeeze(S(i, :, 3));
        J((i + 1):(i + 4), :, :) = 0;
      end
      plot((RI + GI + BI) / 3, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1, 'LineStyle', '-.'); hold on;
      plot(RI, 'Color', 'r', 'LineWidth', 1); hold on;
      plot(GI, 'Color', [0, 0.5, 0], 'LineWidth', 1);
      plot(BI, 'Color', 'b', 'LineWidth', 1);
      plot((-1) * (double(GI) - double(BI)), 'Color', [0, 0.5, 0.5], 'LineWidth', 1);
      plot((-1) * double(RI - GI - BI), 'Color', 'k', 'LineWidth', 1); 
      box off; xlim([1, 2 * nImageHalfWidth]); ylim([-96, 256]);
      if j == 1
        title(sprintf('day %d', nDateDif), 'FontWeight', 'normal');
      end
      j = j + 1;
    end
    subplot(4, 4, 1);
    imshow(J);
    subplot(4, 4, 2);
    X = abs(H(:, :, 2) - H(:, :, 3)) < 10;
    imshow(X);
    subplot(4, 4, 3);
    % X = H(:, :, 1) - H(:, :, 2) - H(:, :, 3) > 0; % dRED1
    X = 2 * H(:, :, 1) - H(:, :, 2) - H(:, :, 3) > 80; % dRED2
    imshow(X);

    H = abs(H(:, :, 2) - H(:, :, 3)) < 10 & 2 * H(:, :, 1) - H(:, :, 2) - H(:, :, 3) > 80;

    % circle limit
    bCircleLimit = 1; % must be 1 always
    if bCircleLimit == 1
      x = sum(H, 1); 
      y = sum(H, 2); 
      x = filtfilt(fb, fa, x);
      y = filtfilt(fb, fa, y);
      % min to max
      % h = 5 / 2;
      % ix0 = find(x > h, 1, 'first');
      % ix1 = find(x > h, 1, 'last');
      % iy0 = find(y > h, 1, 'first');
      % iy1 = find(y > h, 1, 'last');
      % ix = (ix1 - ix0) / 2 + ix0;
      % iy = (iy1 - iy0) / 2 + iy0;
      % max
      % [~, ix] = max(x((nImageHalfWidth - 50):(nImageHalfWidth + 50)));
      % [~, iy] = max(y((nImageHalfWidth - 50):(nImageHalfWidth + 50)));
      % ix = ix + nImageHalfWidth - 50;
      % iy = iy + nImageHalfWidth - 50;
      % center
      % ix = nImageHalfWidth;
      % iy = nImageHalfWidth;

      % cx = ix - nImageHalfWidth;
      % cy = iy - nImageHalfWidth;

      % assume ulcer in the center
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
  
      pR = 5:5:nImageHalfWidth;
      nR = length(pR);
      S = zeros(nR, 1);
      for iR = 1:nR
        R = pR(iR);
        s = sqrt(((-nImageHalfWidth:nImageHalfWidth) - cx) .^ 2 + ((-nImageHalfWidth:nImageHalfWidth)' - cy) .^ 2) < R;
        s = H .* s;
        S(iR) = sum(s(:));
      end
      xR = 10;
      dS = [0; diff(S)];
      hS = dS > xR;
      i = find(hS > 0, 1, 'first');
      if i < 20
        iR = find(hS(i:end) < 1, 1, 'first') + i;
      else
        iR = 1;
      end
      R = pR(iR);
  
      M = sqrt(((-nImageHalfWidth:nImageHalfWidth) - cx) .^ 2 + ((-nImageHalfWidth:nImageHalfWidth)' - cy) .^ 2) < R;
      C = sqrt(((-nImageHalfWidth:nImageHalfWidth) - cx) .^ 2 + ((-nImageHalfWidth:nImageHalfWidth)' - cy) .^ 2) > (R - 5);

      H = M .* H;
    end

    subplot(4, 4, 4);
    imshow(H + C & M);
    title(sprintf('%1.4f', sum(H(:)) / (nImageHalfWidth .^ 2)), 'FontWeight', 'normal');
  
    % save image
    aDir = support_fname({aPath, 'leprosy', '_analysis', 'slicing_plots', aSubject});
    if ~exist(aDir, 'dir')
      mkdir(aDir);
    end
    aFilename = support_fname({aDir, aFile});
    print(hFigure, aFilename, '-dpng', '-r300');
    close(hFigure);  
  end
  o = 0;
end

end % end

%-------------------------------------------------------------------------------




