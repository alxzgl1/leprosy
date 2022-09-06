%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_ulcer_detection()

clc;

% segmentatyion parameters
threshold = 0.08;

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
  for iFile = 1:12 % 1:nFiles
    aFile = tFiles{iFile};
    aFilename = support_fname({aSubpath, aSubject, aFile});
    % load image
    I = imread(aFilename);
    nWidth = size(I, 2);
    nHeight = size(I, 1);

    % focus on image center
    d = 250; % in pixels 
    x = nWidth / 2;
    y = nHeight / 2;
    J = I(:, :, :); 
    J = J((y - d):(y + d), (x - d):(x + d), :);

    subplot(4, 6, 1 + (iFile - 1) * 2);
    imshow(J);

    s = J(:, :, 2);
    s = double(s(:));
    bins = (0:255)';
    H = histc(s, bins);

    subplot(4, 6, 2 + (iFile - 1) * 2);
    bar(bins, H); hold on;
    axis square;


    % plot
    % imshow(I);
    % title('Original image', 'FontWeight', 'normal');
    % hold on;
    % line([nWidth / 2, nWidth / 2], [1, nHeight]); % center line
    % line([1, nWidth], [nHeight / 2, nHeight / 2]); % center line


  end

  return

end


end % end

%-------------------------------------------------------------------------------




