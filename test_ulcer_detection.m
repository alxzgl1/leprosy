%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_ulcer_detection()

clc;

% parameters
bCutImage = 1;
nImageHalfWidth = 350; % in pixels

nSuperPixelClusters = 7; % 500 (uncut image)

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
		
		% cut image
		if bCutImage == 1
      d = nImageHalfWidth;
			x = nWidth / 2;
			y = nHeight / 2;
			I = I(:, :, :); 
			I = I((y - d):(y + d), (x - d):(x + d), :);
		end

    % superpixel
    bSuperPixel = 1;
    if bSuperPixel == 1
      [L, N] = superpixels(I, nSuperPixelClusters);
      BW = boundarymask(L);
      subplot(2, 3, 1);
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
      
      subplot(2, 3, 2);
      imshow(outputImage, 'InitialMagnification', 67);

      % subtract "red - green - blue"
      G = outputImage(:, :, 1) - outputImage(:, :, 2) - outputImage(:, :, 3);
      Q = double(G);

      % get maximum and spot it
      bSpotMaximum = 0;
      if bSpotMaximum == 1
        [~, iMax] = max(G(:));
        x = zeros(size(G(:)));
        x(iMax) = 128;
        y = reshape(x, size(G, 1), size(G, 2));
        [iRow, iCol] = find(y == 128);
        outputImage((iRow - 10):(iRow + 10), (iCol - 10):(iCol + 10), :) = 255;
      end

      B = uint8(Q ~= 0);
      outputImage = outputImage .* repmat(B, 1, 1, 3);
      
      subplot(2, 3, 5);
      imshow(outputImage, 'InitialMagnification', 67); 

      subplot(2, 3, 3);
      imshow(G, 'InitialMagnification', 67); 

      subplot(2, 3, 4);
      imshow(Q, 'InitialMagnification', 67); 

      return
    end
  end
 
  return
end


end % end

%-------------------------------------------------------------------------------




