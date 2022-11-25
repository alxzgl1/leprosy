%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_RGB_decomposition()

% Step 1: Read in the Color Image and Convert it to Grayscale
% rgb = imread('pears.png');

I = imread('d:\data\leprosy\TABLE_Aranz_Image\S-01\S-1_20201002.png');

% cut image
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

D = [8, 8];

J_R = medfilt2(I(:, :, 1), D);
J_G = medfilt2(I(:, :, 2), D);
J_B = medfilt2(I(:, :, 3), D);

J = cat(3, J_R, J_G, J_B);

H = abs(J(:, :, 2) - J(:, :, 3)) < 10 & J(:, :, 1) - J(:, :, 2) - J(:, :, 3) > 0;

iSlice = 350;
I((iSlice + 1):(iSlice + 3), :, :) = 0;
J((iSlice + 1):(iSlice + 3), :, :) = 0;

figure;

subplot(2, 3, 1); imagesc(I);
subplot(2, 3, 2); imagesc(J);
subplot(2, 3, 3); imagesc(H);

subplot(2, 3, 4); plot(I(iSlice, :, 1), 'r'); hold on; plot(I(iSlice, :, 2), 'g'); plot(I(iSlice, :, 3), 'b'); 
xlim([1, 700]); ylim([0, 255]); box off;

subplot(2, 3, 5); plot(J(iSlice, :, 1), 'r'); hold on; plot(J(iSlice, :, 2), 'g'); plot(J(iSlice, :, 3), 'b'); 
xlim([1, 700]); ylim([0, 255]); box off;

subplot(2, 3, 6); plot(J(iSlice, :, 2) - J(iSlice, :, 3), 'k'); 
xlim([1, 700]); ylim([0, 255]); box off;

end % end

%-------------------------------------------------------------------------------