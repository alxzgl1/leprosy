%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_k_means()

% https://uk.mathworks.com/help/images/color-based-segmentation-using-k-means-clustering.html

% parameters
bCutImage = 1;
nImageHalfWidth = 350; % in pixels | 350 (default) | stable parameter

% get path
aPath = support_get_path();
aSubpath = support_fname({aPath, 'leprosy', 'TABLE_Aranz_Image'});

% get subjects
tSubjects = support_get_subjects(aSubpath);
nSubjects = length(tSubjects);

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

    % subtraction
    % bSubtraction = 0;
    % if bSubtraction == 1
    %   x = J(:, :, 1) - I(:, :, 2) - I(:, :, 3);
    %   J = cat(3, x, 0 * x, 0 * x);
    % end

     % open figure
    hFigure = figure; 
    set(hFigure, 'NumberTitle', 'off', 'Position', [0, 0, 1920, 1080] / 2.0, 'MenuBar', 'none', 'Resize', 'off', 'Visible', 'off'); 

    % plot
    subplot(2, 4, 1); imshow(I);
    subplot(2, 4, 2); imshow(J);

    he = J;
    
    numColors = 3;
    L = imsegkmeans(he, numColors);
    B = labeloverlay(he, L);
    subplot(2, 4, 3); imshow(B);
    title("Labeled Image RGB");
    
    lab_he = rgb2lab(he);
    
    ab = lab_he(:,:,2:3);
    ab = im2single(ab);
    pixel_labels = imsegkmeans(ab,numColors,NumAttempts=3);
    
    B2 = labeloverlay(he,pixel_labels);
    subplot(2, 4, 4); imshow(B2);
    title("Labeled Image a*b*");
    
    
    mask1 = pixel_labels == 1;
    cluster1 = he.*uint8(mask1);
    subplot(2, 4, 5); imshow(cluster1);
    title("Objects in Cluster 1");
    
    
    mask2 = pixel_labels == 2;
    cluster2 = he.*uint8(mask2);
    subplot(2, 4, 6); imshow(cluster2);
    title("Objects in Cluster 2");
    
    
    mask3 = pixel_labels == 3;
    cluster3 = he.*uint8(mask3);
    subplot(2, 4, 7); imshow(cluster3)
    title("Objects in Cluster 3");
    
    
    L = lab_he(:,:,1);
    L_blue = L.*double(mask1);
    L_blue = rescale(L_blue);
    idx_light_blue = imbinarize(nonzeros(L_blue));
    
    blue_idx = find(mask1);
    mask_dark_blue = mask1;
    mask_dark_blue(blue_idx(idx_light_blue)) = 0;
    
    blue_nuclei = he.*uint8(mask_dark_blue);
    subplot(2, 4, 8); imshow(blue_nuclei)
    title("Blue Nuclei")

    % save image
    aDir = support_fname({aPath, 'leprosy', '_analysis', 'k_means_plots', aSubject});
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