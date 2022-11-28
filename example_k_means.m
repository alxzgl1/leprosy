%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function example_k_means()

% https://uk.mathworks.com/help/images/color-based-segmentation-using-k-means-clustering.html

% Step 1: Read in the Color Image and Convert it to Grayscale
% rgb = imread('pears.png');

rgb = imread('d:\data\leprosy\TABLE_Aranz_Image\S-01\S-1_20200921.png');

% cut image
I = rgb;
nWidth = size(I, 2);
nHeight = size(I, 1);
nImageHalfWidth = 350;
% cut image
bCutImage = 1;
if bCutImage == 1
  d = nImageHalfWidth;
	x = nWidth / 2;
	y = nHeight / 2;
	I = I(:, :, :); 
	I = I((y - d):(y + d), (x - d):(x + d), :);
end
subplot(2, 4, 1); imshow(I);

% median filter
D = [32, 32];
J_R = medfilt2(I(:, :, 1), D);
J_G = medfilt2(I(:, :, 2), D);
J_B = medfilt2(I(:, :, 3), D);
J = cat(3, J_R, J_G, J_B);

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

end % end

%-------------------------------------------------------------------------------