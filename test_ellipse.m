%-------------------------------------------------------------------------------
% Function
% LINK:
%   * https://uk.mathworks.com/discovery/image-segmentation.html
%   * https://uk.mathworks.com/help/vision/ug/segment-3d-brain-tumor-using-deep-learning.html
%-------------------------------------------------------------------------------
function test_ellipse()

clc;

nWidth = 100;
nHeight = 100; 
I = zeros(nHeight, nWidth);

% nCenter = [ix, iy];
% nRadius = [rx, ry];

eccentricity = 1;
% 
% x1 = 1;
% x2 = 10;
% y1 = 1;
% y2 = 10;
% 
% 
% 
%     a = (1/2) * sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
%     b = a * sqrt(1 - eccentricity ^ 2);
%     t = linspace(0, 2 * pi, 100); % Absolute angle parameter
%     X = a * cos(t);
%     Y = b * sin(t);
%     % Compute angles relative to (x1, y1).
%     angles = atan2(y2 - y1, x2 - x1);
%     x = (x1 + x2) / 2 + X * cos(angles) - Y * sin(angles);
%     y = (y1 + y2) / 2 + X * sin(angles) + Y * cos(angles);
% 
%     plot(x, y, '.');

% OLD
bOLD = 1;
if bOLD == 1
  phi = -45; 
  rx = 1;
  ry = 2;
  t = 0:0.1:2*pi;
  phi_rad = (phi / 180) * pi;
  x = rx * cos(t) * cos(phi_rad) - ry * sin(t) * sin(phi_rad);
  y = rx * cos(t) * sin(phi_rad) + ry * sin(t) * cos(phi_rad);
  plot(x, y, '.');
end


% for i = 
plot(x, y, 'Color', 'r', 'LineWidth', 3);


end % end

%-------------------------------------------------------------------------------




