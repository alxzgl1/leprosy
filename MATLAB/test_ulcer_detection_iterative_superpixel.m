%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_ulcer_detection_iterative_superpixel()

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

      % iterate clusters
      nClusters = 25;
      nOptimalClusters = 0;
      for iCluster = 1:nClusters
        % do clutering
        [L, N] = superpixels(I, iCluster);
        BW = boundarymask(L);
        % init image
        outputImage = zeros(size(I), 'like', I);
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
        % do colour difference
        G = outputImage(:, :, 1) - outputImage(:, :, 2) - outputImage(:, :, 3);
        Q = double(G);
        % quit when cluster is found
        if sum(Q(:)) > 0
          nOptimalClusters = iCluster;
          break
        end
      end

      % * panel A *

      % redo clustering
      [L, N] = superpixels(I, nOptimalClusters);
      BW = boundarymask(L);
      % init image
      outputImage = zeros(size(I), 'like', I);
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
      % color difference
      G = outputImage(:, :, 1) - outputImage(:, :, 2) - outputImage(:, :, 3);
      G = G > 0;
      G = uint8(repmat(G, 1, 1, 3) * 255);
      % show image
      subplot(2, 4, 1);
      imshow(imoverlay(I, BW, 'cyan'), 'InitialMagnification', 67);
      subplot(2, 4, 2);
      imshow(outputImage, 'InitialMagnification', 67);
      subplot(2, 4, 3);
      imshow(G);
      GI_A = uint8(G > 0) .* outputImage;
      subplot(2, 4, 4);
      imshow(GI_A);

      % * panel B *

      % redo clustering
      nMaxClusters = 65;
      [L, N] = superpixels(I, nMaxClusters);
      BW = boundarymask(L);
      % init image
      outputImage = zeros(size(I), 'like', I);
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
      % show image
      subplot(2, 4, 5);
      imshow(imoverlay(I, BW, 'cyan'), 'InitialMagnification', 67);
      subplot(2, 4, 6);
      imshow(outputImage, 'InitialMagnification', 67);
      subplot(2, 4, 7);
      G = outputImage(:, :, 1) - outputImage(:, :, 2) - outputImage(:, :, 3);
      U = unique(G(:));
      G = G > U(3); % > 0 (default) <- tunning parameter
      G = uint8(repmat(G, 1, 1, 3) * 255);
      imshow(G);
      GI_B = uint8(G > 0) .* outputImage;
      subplot(2, 4, 8);
      imshow(GI_B);

      return
    end
  end
 
  return
end


end % end

%-------------------------------------------------------------------------------




