%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function example_circular_objects_detection()

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
rgb = I;

[centers, radii] = imfindcircles(rgb, [50 200], 'ObjectPolarity', 'dark', ...
    'Sensitivity', 0.95);

radii

imshow(rgb)
h = viscircles(centers,radii);

% https://uk.mathworks.com/help/images/detect-and-measure-circular-objects-in-an-image.html

end % end

%-------------------------------------------------------------------------------