%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_ulcer_slicing_plots_w_filter()

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
[fb, fa] = butter(4, 0.05, 'low');

% loop subjects
for iSubject = 1:nSubjects
  aSubject = tSubjects{iSubject};
  % status
  fprintf(1, '%s\n', aSubject);

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
  nDateRef = [];
  nFiles = length(tFiles);
  for iFile = 1:nFiles
    aFile = tFiles{iFile};
    aFilename = support_fname({aSubpath, aSubject, aFile});

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
    bMedianFilter = 0;
    if bMedianFilter == 1
      D = [20, 20];
      J_R = medfilt2(I(:, :, 1), D);
      J_G = medfilt2(I(:, :, 2), D);
      J_B = medfilt2(I(:, :, 3), D);
      J = cat(3, J_R, J_G, J_B);
      H = double(J);
      S = J;
    end

    % filter
    bLP_Filter = 1;
    if bLP_Filter == 1
      J_R = I(:, :, 1);
      J_G = I(:, :, 2);
      J_B = I(:, :, 3);
      % smooth rows
      J_R = uint8(filtfilt(fb, fa, double(J_R)')');
      J_G = uint8(filtfilt(fb, fa, double(J_G)')');
      J_B = uint8(filtfilt(fb, fa, double(J_B)')');
      % smooth cols
      J_R = uint8(filtfilt(fb, fa, double(J_R)));
      J_G = uint8(filtfilt(fb, fa, double(J_G)));
      J_B = uint8(filtfilt(fb, fa, double(J_B)));
      % init
      J = cat(3, J_R, J_G, J_B);
      H = double(J);
      S = double(J);
    end

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
    x = abs(H(:, :, 2) - H(:, :, 3)) < 10;
    imshow(x);
    subplot(4, 4, 3);
    x = H(:, :, 1) - H(:, :, 2) - H(:, :, 3) > 0;
    imshow(x);
    subplot(4, 4, 4);
    H = abs(H(:, :, 2) - H(:, :, 3)) < 10 & H(:, :, 1) - H(:, :, 2) - H(:, :, 3) > 0;
    imshow(H);
  
    % save image
    aDir = support_fname({aPath, 'leprosy', '_analysis', 'slicing_plots_w_filter', aSubject});
    if ~exist(aDir, 'dir')
      mkdir(aDir);
    end
    aFilename = support_fname({aDir, aFile});
    print(hFigure, aFilename, '-dpng', '-r300');
    close(hFigure);  
  end
end

end % end

%-------------------------------------------------------------------------------




