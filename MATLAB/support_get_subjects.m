%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function tSubjects = support_get_subjects(aPath)

tDir = dir(aPath);

tSubjects = [];

k = 1;
for i = 1:length(tDir)
  a = tDir(i).name;
  if contains(a, 'S-')
    tSubjects{k} = a;
    k = k + 1;
  end
end

end % end

%-------------------------------------------------------------------------------