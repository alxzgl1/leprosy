%-------------------------------------------------------------------------------
% Function
% LINK:
%   * https://uk.mathworks.com/discovery/image-segmentation.html
%   * https://uk.mathworks.com/help/vision/ug/segment-3d-brain-tumor-using-deep-learning.html
%-------------------------------------------------------------------------------
function test_ellipse()

clc;

nWidth = 800;
nHeight = 600; 

% ellipse
phi = 135;
rx = 200;
ry = 100;
t = linspace(0, 2 * pi, 1000); 
phi_rad = (phi / 180) * pi;
x = rx * cos(t) * cos(phi_rad) - ry * sin(t) * sin(phi_rad);
y = rx * cos(t) * sin(phi_rad) + ry * sin(t) * cos(phi_rad);

% plot
subplot(2, 2, 1);
plot(x, y, '.'); box off;
axis square;
xlim([-max(rx, ry), max(rx, ry)]); 
ylim([-max(rx, ry), max(rx, ry)]);

% make image circle
CI = zeros(nHeight, nWidth);

gx = round(x + nWidth / 2);
gy = round(y + nHeight / 2);

for i = 1:length(t)
  CI(gy(i), ((gx(i) - 1):(gx(i) + 1))) = 1;
end

gy_min = min(gy);
gy_max = max(gy);
for i = gy_min:gy_max
  K = find(CI(i, :) > 0);
  if ~isempty(K)
    CI(i, K(1):K(end)) = 1;
  end
end

% plot
subplot(2, 2, 2);
imagesc(CI); 
axis square;
set(gca, 'YDir', 'normal');

end % end

%-------------------------------------------------------------------------------




