%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function test_watershed_segmentation()

% Step 1: Read in the Color Image and Convert it to Grayscale
% rgb = imread('pears.png');

rgb = imread('c:\Users\alexa\Downloads\S-1_20200921.png');

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
rgb = I;

I = rgb2gray(rgb);
subplot(2, 4, 1); imshow(I)
text(732,501,'Image courtesy of Corel(R)',...
     'FontSize',7,'HorizontalAlignment','right')

% Step 2: Use the Gradient Magnitude as the Segmentation Function
gmag = imgradient(I);
subplot(2, 4, 2); imshow(gmag,[])
title('Gradient Magnitude')

% Step 3: Mark the Foreground Objects
se = strel('disk',20);
Io = imopen(I,se);
subplot(2, 4, 3); imshow(Io)
title('Opening')

Ie = imerode(I,se);
Iobr = imreconstruct(Ie,I);
subplot(2, 4, 4); imshow(Iobr)
title('Opening-by-Reconstruction')

Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
subplot(2, 4, 5); imshow(Iobrcbr)
title('Opening-Closing by Reconstruction')

fgm = imregionalmax(Iobrcbr);
subplot(2, 4, 6); imshow(fgm)
title('Regional Maxima of Opening-Closing by Reconstruction')

% Step 4: Compute Background Markers
bw = imbinarize(Iobrcbr);
subplot(2, 4, 7); imshow(bw)
title('Thresholded Opening-Closing by Reconstruction')

end % end

%-------------------------------------------------------------------------------