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

% smoothing filter
[fb, fa] = butter(2, 0.1, 'low');

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

    % superpixel
    bSuperPixel = 1;
    if bSuperPixel == 1
      [L, N] = superpixels(I, 500);
      BW = boundarymask(L);
      subplot(2, 2, 1);
      imshow(imoverlay(I, BW, 'cyan'), 'InitialMagnification', 67);

      % https://uk.mathworks.com/help/images/ref/superpixels.html

      outputImage = zeros(size(I),'like',I);
      idx = label2idx(L);
      numRows = size(I,1);
      numCols = size(I,2);
      for labelVal = 1:N
        redIdx = idx{labelVal};
        greenIdx = idx{labelVal} + numRows * numCols;
        blueIdx = idx{labelVal} + 2 * numRows * numCols;
        outputImage(redIdx) = mean(I(redIdx));
        outputImage(greenIdx) = mean(I(greenIdx));
        outputImage(blueIdx) = mean(I(blueIdx));
      end    
      
      subplot(2, 2, 2);
      imshow(outputImage, 'InitialMagnification', 67);

      % histogram
      subplot(2, 2, 3);
      s = outputImage(:, :, 1);
      s = double(s(:));

      p = sort(s, 'descend');
      p = unique(p);

      % ---- HERE ---- need to be fixed

      Y = uint8(zeros(size(outputImage)));
      for i = 1:8
        Y(outputImage == p(i)) = 1;
      end
      imshow(imoverlay(Y, BW, 'cyan'), 'InitialMagnification', 67);

%       bins = (0:255)';
%       H = histc(s, bins) / length(bins);
%       % H = filtfilt(fb, fa, H);
%       plot(bins, H); hold on; box off;
%       xlim([0, 255]);
%       axis square;


      return
    end

    % focus on image center
    d = 250; % in pixels 
    x = nWidth / 2;
    y = nHeight / 2;
    J = I(:, :, :); 
    J = J((y - d):(y + d), (x - d):(x + d), :);

    subplot(4, 6, 1 + (iFile - 1) * 2);

    % remove RBG
    J(:, :, 1) = uint8(J(:, :, 1) > 100) .* J(:, :, 1);
    J(:, :, 2) = uint8(J(:, :, 2) > 100) .* J(:, :, 2);

    imshow(J);

    subplot(4, 6, 2 + (iFile - 1) * 2);

    % loop RGB
    tColors = {'r', 'g', 'b'};
    for i = 1:3
      s = J(:, :, i);
      s = double(s(:));
      bins = (0:255)';
      H = histc(s, bins) / length(bins);
      % H = filtfilt(fb, fa, H);
      plot(bins, H, 'Color', tColors{i}); hold on; box off;
      xlim([0, 255]);
      axis square;
    end

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




