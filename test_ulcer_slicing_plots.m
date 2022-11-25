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
% [fb, fa] = butter(4, 0.3, 'low');

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
    D = [32, 32];
    J_R = medfilt2(I(:, :, 1), D);
    J_G = medfilt2(I(:, :, 2), D);
    J_B = medfilt2(I(:, :, 3), D);
    J = cat(3, J_R, J_G, J_B);
    H = J;

    % open figure
    hFigure = figure; 
    set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 2.0, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 

    % decompose
    j = 1;
    for i = 50:50:(2 * nImageHalfWidth)
      subplot(4, 4, j + 2);
      bCols = 0;
      if bCols == 1
        RI = squeeze(I(:, i, 1));
        GI = squeeze(I(:, i, 2));
        BI = squeeze(I(:, i, 3));
        J(:, (i + 1):(i + 4), :) = 0;
      else
        RI = squeeze(I(i, :, 1));
        GI = squeeze(I(i, :, 2));
        BI = squeeze(I(i, :, 3));
        J((i + 1):(i + 4), :, :) = 0;
      end
      plot((double(RI) + double(GI) + double(BI)) / 3, 'Color', [0.5, 0.5, 0.5]); hold on;
      plot(RI, 'Color', 'r'); 
      plot(GI, 'Color', [0, 0.5, 0]);
      plot(BI, 'Color', 'b');
      plot((-1) * (double(GI) - double(BI)), 'Color', [0, 0.5, 0.5]);
      plot((-1) * double(RI - GI - BI), 'Color', 'k'); 
      box off; xlim([1, 2 * nImageHalfWidth]); ylim([-96, 256]);
      if j == 1
        title(sprintf('day %d', nDateDif), 'FontWeight', 'normal');
      end
      j = j + 1;
    end
    subplot(4, 4, 1);
    imshow(J);
    subplot(4, 4, 2);
    H = abs(H(:, :, 2) - H(:, :, 3)) < 10 & H(:, :, 1) - H(:, :, 2) - H(:, :, 3) > 0;
    imshow(H);
  
    % save image
    aDir = support_fname({aPath, 'leprosy', '_analysis', 'slicing_plots', aSubject});
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




