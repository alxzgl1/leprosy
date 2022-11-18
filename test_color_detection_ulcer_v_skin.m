%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_color_detection_ulcer_v_skin()

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
[fb, fa] = butter(4, 0.3, 'low');

% loop subjects
for iSubject = 1:1 % nSubjects
  % aSubject = tSubjects{iSubject};

  aSubject = tSubjects{contains(tSubjects, 'S-175')};

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

  % decompose
  bDecompose = 1;
  if bDecompose == 1
    j = 1;
    J = I;
    for i = 1:25:700
      subplot(4, 7, j);
      bCols = 1;
      if bCols == 1
        RI = squeeze(I(:, i, 1));
        GI = squeeze(I(:, i, 2));
        BI = squeeze(I(:, i, 3));
        J(:, i, :) = 0;
      else
        RI = squeeze(I(i, :, 1));
        GI = squeeze(I(i, :, 2));
        BI = squeeze(I(i, :, 3));
        J(i, :, :) = 0;
      end
      bGB_diff = 0;
      if bGB_diff == 1
        plot(RI, 'Color', 'r'); hold on;
        plot(GI + BI, 'Color', [0, 0.5, 0.5]);
        plot(RI - GI - BI, 'Color', 'k');
      else
        plot(RI, 'Color', 'r'); hold on;
        plot(GI, 'Color', [0, 0.5, 0]);
        plot(BI, 'Color', 'b');
        plot(RI - GI - BI, 'Color', 'k'); 
      end
      j = j + 1;
    end
    figure;
    imshow(J);
    figure;
    imshow(I(:, :, 1) - I(:, :, 2) - I(:, :, 3));
  end

end

end % end

%-------------------------------------------------------------------------------




