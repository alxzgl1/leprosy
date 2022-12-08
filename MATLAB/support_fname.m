%-------------------------------------------------------------------------------
% Function
%-------------------------------------------------------------------------------
function aFilename = support_fname(tArgs)

aFilename = tArgs{1};
for nArg = 2:length(tArgs)
	aFilename = sprintf('%s%s%s', aFilename, filesep, tArgs{nArg});
end

end % end 

%-------------------------------------------------------------------------------