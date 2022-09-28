%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_simple_color_difference()

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
[fb, fa] = butter(2, 0.1, 'low');

% loop subjects
for iSubject = 5:nSubjects
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
    G = I(:, :, 1) - I(:, :, 2) - I(:, :, 3);
    G = repmat(G, 1, 1, 3);

    GF = I(:, :, 1) - I(:, :, 2) - I(:, :, 3);
    GF = medfilt2(GF, [8, 8]);

    h = 32;
    Q = uint8((GF > h) * 255);

    % show image
    subplot(2, 4, 1);
    imshow(I);
    subplot(2, 4, 2);
    imshow(G);
    subplot(2, 4, 3);
    imshow(Q);
    QI = uint8(Q > 0) .* I;
    subplot(2, 4, 4);
    imshow(QI);

    return

  end
 
  return
end


end % end

%-------------------------------------------------------------------------------




