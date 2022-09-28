%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_plot_images_over_time()

clc;

bHistogram = 0;

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

    % focus on image center
    d = 250; % in pixels 
    x = nWidth / 2;
    y = nHeight / 2;
    J = I(:, :, :); 
    J = J((y - d):(y + d), (x - d):(x + d), :);

    % histogram
    if bHistogram == 1
      % image
      subplot(4, 6, 1 + (iFile - 1) * 2);
      % remove RBG
      % J(:, :, 1) = uint8(J(:, :, 1) > 100) .* J(:, :, 1);
      % J(:, :, 2) = uint8(J(:, :, 2) > 100) .* J(:, :, 2);
      imshow(J);

      % histogram
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
    else

      % image
      subplot(3, 4, iFile);
      % remove RBG
      % J(:, :, 1) = uint8(J(:, :, 1) > 100) .* J(:, :, 1);
      % J(:, :, 2) = uint8(J(:, :, 2) > 100) .* J(:, :, 2);
      imshow(J);
      title(sprintf('week #%d', iFile), 'FontSize', 10, 'FontWeight', 'normal');
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




