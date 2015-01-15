%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (fakeMat)
%
% This function can be used to generate matrices readable by OpenBR.
% Use this function to test OpenBR's output and assumptions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % These are example values of the function's inputs
% filename = 'check';
% mtx = %[-Inf -Inf -Inf ; 10 10 5 ; 10 5 10];
% mask = [-1 127 127 ; 127 -1 127 ; 127 127 -1];
% dimensions = [3 3];

%
% NOTE: A mask value of -1 indicates Genuine (match), while 127 indicates
% Impostor (non-match). Values of 0 are ignored.
%

function [] = fakeMat(filename, mtx, mask, dimensions)

% Open matrix file
mtx_file = fopen([filename '.mtx'], 'w');
% Write header with appropriate matrix dimensions
fwrite(mtx_file, sprintf(['S2\n\n\nMF ' num2str(dimensions(1)) ' ' num2str(dimensions(2)) '\n']));
% Write data as single precision float
fwrite(mtx_file, mtx, 'single');
% Close the file
fclose(mtx_file);

% Open matrix file
mask_file = fopen([filename '.mask'], 'w');
% Write header with appropriate matrix dimensions
fwrite(mask_file, sprintf(['S2\n\n\nMB ' num2str(dimensions(1)) ' ' num2str(dimensions(2)) '\n']));
% Write data as 8-bit integer
fwrite(mask_file, mask, 'int8');
% Close the file
fclose(mask_file);

end