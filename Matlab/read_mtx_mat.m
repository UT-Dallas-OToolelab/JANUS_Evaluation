function mtx_value = read_mtx_mat(mtx_name, n_gal, n_probe)
% Format: Output = read_mtx_mat('name of matrix', # items in gallery,
% # of items in probe)
%====
%Used to read binary matrices for JANUS project. cd to folder containing
%matrix to be analyzed.
%Number of gallery and probe items are obtained by opening .mtx (in
%benchmarks dir) and looking at the numerical values in line 4. The first
%value is the number of probe items. The second value is the number of
%gallery items.

%% Read matrix

%open file
fid = fopen(mtx_name);
%read in parameters
[~,~,machinefmt,~] = fopen(fid);
%Position to end of file, without moving any bytes from the origin
fseek(fid, 0, 'eof');
%get current position in file (will correspond to file length)
filelength = ftell(fid);
%Position fo begining of file, withouth moving any bytes from the origin
fseek(fid, 0, 'bof');


%% Format matrix

%find number of items in the header
header_n = filelength-(n_gal*n_probe*4);
%read in header: get data from file, size of header, precision, skipping no
%bytes, in current machine's format
header = fread(fid,header_n,'uint8=>char',0,machinefmt);
%read in matrix: get data from file, size of data, precision,in current
%machine's format
mtx_value = fread(fid,n_gal*n_probe*4,'float=>double',0,machinefmt);

%% Plot values

%Make histogram trimming off max and min values from the plot, with X number of bins
%%histogram(mtx_value(mtx_value > min(mtx_value) & mtx_value < max(mtx_value)),100);

%Make histogram no trimming
%histogram(mtx_value,100);

%% Close file
fclose(fid);


end
