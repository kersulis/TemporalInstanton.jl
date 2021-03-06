%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\Jonas\Box Sync\Research\LANL Summer 2014\RTS-96
%    Data\Table-01.txt
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

%% Initialize variables.
filename = 'C:\Users\Jonas\Box Sync\Research\LANL Summer 2014\RTS-96 Data\Table-01.txt';
startRow = 5;
endRow = 77;

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%3s%9s%5s%5s%5s%3s%5s%3s%6s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,3,4,5,6,7,8,9,10]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [1,3,4,5,6,7,8,9,10]);
rawCellColumns = raw(:, 2);

%% Allocate imported array to column variable names
BusNo = cell2mat(rawNumericColumns(:, 1));
BusName = rawCellColumns(:, 1);
BusType = cell2mat(rawNumericColumns(:, 2));
LoadMW = cell2mat(rawNumericColumns(:, 3));
LoadMVar = cell2mat(rawNumericColumns(:, 4));
GL = cell2mat(rawNumericColumns(:, 5));
BL = cell2mat(rawNumericColumns(:, 6));
SubArea = cell2mat(rawNumericColumns(:, 7));
Base_kV = cell2mat(rawNumericColumns(:, 8));
Zone = cell2mat(rawNumericColumns(:, 9));

%% Clear temporary variables
clearvars filename startRow endRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns;

%% Begin MATPOWER Data Formatting
% Format bus data for 73 buses
% BusType(:) = 1; % All P and Q are specified, except at...
% BusType(13) = 3; % ...Arne, which is the slack bus.

Vm(1:73,1) = 0; % Voltages are not specified at PQ or slack buses
Va(1:73,1) = 0;

Vmax(1:73,1) = 1.05; % Set voltage limits to typical values
Vmin(1:73,1) = 0.95;

% The following BusData array may be pasted into a MATPOWER case
% definition:
format ShortG % Change format to remove trailing zeros
BusData = [BusNo BusType LoadMW LoadMVar GL BL SubArea Vm Va Base_kV Zone Vmax Vmin]
format
