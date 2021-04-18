function [ plasma ] = get_CHARLIE_nonuniform_plasma( in_type)
% get ne/Te of spatially nonuniform plasma of CHARLIE, accroding to the
% nonuniform ne/Te distribution of CHARLIE(1MHz, 1Pa, 520W). Details are in
% fit_CHARLIE_nonuniform_plasma.m

%% io
% in_type: str
% meaning of characters
% '': output all
% first character: 'z' for axially, or 'r' for radially
% end character: 'a' for average
% number b: position, at b mm
% 'b-c': region from b mm to c mm. must end with 'a'.
% 'pn': equally divided into n part. must end with 'a'.

% plasma: struct
% plasma.ne/Te: norm ne/Te distribution related to the input type
% plasma

%% nonuniform distribution data
ref.ne_r;


%% operation according to input type
if isempty(in_type)
    plasma=ref;
end
type = strsplit(in_type, '_');


switch type
    case 'r0'
        
    otherwise
        error('No such type!')
end


end