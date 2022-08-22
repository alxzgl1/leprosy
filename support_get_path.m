%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function aPath = support_get_path(argin)

% check os
if contains(getenv('OS'), 'Win') 
	aPath = 'z:\zhigaloa\data';
else 
  aPath = '/rds/projects/2018/jenseno-bbfreqtagging/zhigaloa/data';
end

% WARING
% fprintf(1, 'MARK BAD CHANNELS: MEG1843 (noisy), MEG2122 (flat)\n');

% home computer (HP)
if strcmp(getenv('COMPUTERNAME'), 'DESKTOP-31B7OD2')
  aPath = 'g:\\data';
  return
end

% home computer (Lenovo)
if strcmp(getenv('COMPUTERNAME'), 'LAPTOP-TFCS9T7P')
  aPath = 'd:\\data';
  return
end

% desktop
if nargin > 0
  aPath = 'c:\\data';
end

end % end

%-------------------------------------------------------------------------------