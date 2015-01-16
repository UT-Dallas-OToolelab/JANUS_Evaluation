function mask_vals = read_mask_mat(mask_name,n_gal,n_probe)
% Format: output = read_mask_mat('name of mask matrix',number items in gallery set, number of items in probe set)  
% ====
% Used to read binary masks for JANUS project: cd to folder that contains mask to be analyzed 
% Number of gallery and probe items are obtained by opening .mask file (in
% benchmarks dir) and looking at the numerical values in line 4. The first
% value is the number of probe times. The second value is the number of
% gallery items.
% ====

%% Read mask

%open file
fid = fopen(mask_name);
%read in parameters
[~,~,machinefmt,~] = fopen(fid);
%Position to end of file, without moving any bytes from the origin
fseek(fid, 0, 'eof');
%get current position in file (will correspond to file length)
filelength = ftell(fid);
%Position fo begining of file, withouth moving any bytes from the origin
fseek(fid, 0, 'bof');

%% Format mask

%find number of items in the header
header_n = filelength-(n_gal*n_probe);
%read in header: get data from file, size of header, precision of number of
%byes, skipping no bytes, in current machine's format
header = fread(fid,header_n,'uint8=>char',0,machinefmt);
%read in mask: get data from file, size of data, precision,skipping no
%bytes, in current machine's format
mask_vals = fread(fid,n_gal*n_probe,'uint8=>double',0,machinefmt);

%% Close file
fclose(fid);

end
