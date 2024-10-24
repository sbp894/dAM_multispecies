function GER122E1 = read_txt_data(filename, dataLines)
%IMPORTFILE Import data from a text file
%  GER122E1 = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  GER122E1 = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  GER122E1 = importfile("D:\Dropbox\Pitts_files\Aravind_ExpAM\Exports\DoDAM1k\GER12-2_E1", [7, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 01-Oct-2022 15:06:08

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [7, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["RUN2LEVEL", "Var2", "Var3"];
opts.SelectedVariableNames = "RUN2LEVEL";
opts.VariableTypes = ["double", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var2", "Var3"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var2", "Var3"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "RUN2LEVEL", "TrimNonNumeric", true);
opts = setvaropts(opts, "RUN2LEVEL", "ThousandsSeparator", ",");

% Import the data
GER122E1 = readtable(filename, opts);
GER122E1 = table2array(GER122E1);
end