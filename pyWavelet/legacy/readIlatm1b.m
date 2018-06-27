function atmPre = readIlatm1b(filePath)
% readIlatm1b - This function loads the atm file given the specified
% file path and file name.
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    filePath - The file path indicating the location of the specified ATM file
%
% Outputs:
%    atmPre - A structure containing the loaded atm parameters.
%
% Example:
% -
%
% Other m-files required: none
%
% Subfunctions: none
%
% MAT-files required: none
%
% See also: getAtmDataForFrame
%
%
% Author: Thomas Newman, Ph.D., Earth Sciences
%
% February 2015; Last revision: Base.

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% BEGIN CODE
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%==========================================================================
%% Check if the data is in binary (.qi) or HDF5 (.h5) form.
%==========================================================================
if ~isempty(strfind(filePath,'.qi'))
    dataType = 'binary';
    
elseif ~isempty(strfind(filePath,'.h5'))
    dataType = 'hdf5';
end

%==========================================================================
%% If the ATM data is in binary format.
%==========================================================================
switch dataType
    case 'binary'
        
        %--------------------------------------------------------------------------
        % Reads in data.
        %--------------------------------------------------------------------------
        fclose('all');
        % Closes all open files.
        
        fid = fopen(filePath,'r');
        % Opens the file, filePath, for binary read access.
        
        numberOfWords = fread(fid,[1,1],'int')/4;
        % The number of words specified at the beginning of the binary file.
        
        if numberOfWords>1000
            fclose(fid);
            fid = fopen(filePath,'r','b');
            numberOfWords = fread(fid,[1,1],'int')/4;
        end
        % Check endian of input data.
        
        blankArray    = fread(fid,[1,numberOfWords-1],'int');
        % Reads the blank array.
        
        initialWord   = fread(fid,[1,1],'int');
        % The initial word.
        
        skipBytes     = fread(fid,[1,1],'int');
        % The numer of bytes to skip.
        
        frewind(fid);
        % Rewinds the pointers to the beginning of the file.
        
        fseek(fid, skipBytes, 'bof');
        % Move to specified position in file.
        
        Data = fread(fid,[12,inf],'int');
        % Reads the data array
        
        fclose(fid);
        % Closes an open file.
        
        %--------------------------------------------------------------------------
        % Assigns data to the atmPre data structure.
        %--------------------------------------------------------------------------
        % Data(1 ,:)/1000;
        % (1) Relative Time (seconds from start of the file)
        
        atmPre.latitude = Data(2 ,:)/1e6;
        % (2) Laser Spot Latitude (decimal degrees) *
        
        atmPre.longitude = Data(3 ,:)/1e6;
        % (3) Laser Spot Longitude (decimal degrees) *
        
        atmPre.elevation = Data(4 ,:)/1000;
        % (4) Elevation (meters) *
        
        % Data(5 ,:);
        % (5) Start Pulse signal Strength (relative integer)
        
        % Data(6 ,:);
        % (6) Reflected Laser Signal Strength (relative integer)
        
        atmPre.scanAzimuth  = Data(7 ,:)/1000;
        % (7) Scan Azimuth (degrees) *
        
        atmPre.pitch  = Data(8 ,:)/1000;
        % (8) Pitch (degrees) *
        
        atmPre.roll  = Data(9 ,:)/1000;
        % (9) Roll (degrees) *
        
        % Data(10,:)/10;
        % (10) GPS PDOP (dilution of precision)
        
        % Data(11,:);
        % (11) Laser received pulse width (digitizer samples)
        
        atmPre.gpsTimeHHMMSS = Data(12,:);
        % (12) GPS time packed,
        % Example: 153320.100 = 15 hours 33 minutes 20 seconds 100 milliseconds.
        
end

%==========================================================================
%% If the ATM data is in HDF5 format.
%==========================================================================
switch dataType
    case 'hdf5'
        
        atmPre.latitude = h5read(filePath,'/latitude').';
        % Laser Spot Latitude (decimal degrees) *
        
        atmPre.longitude = h5read(filePath,'/longitude').';
        % Laser Spot Longitude (decimal degrees) *
        
        atmPre.elevation = h5read(filePath,'/elevation').';
        % Elevation (meters) *
        
        atmPre.scanAzimuth = h5read(filePath,'/instrument_parameters/azimuth').';
        % Scan Azimuth (degrees) *
        
        atmPre.pitch = h5read(filePath,'/instrument_parameters/pitch').';
        % Pitch (degrees) *
        
        atmPre.roll = h5read(filePath,'/instrument_parameters/roll').';
        % Roll (degrees) *
        
        atmPre.gpsTimeHHMMSS = 1000*h5read(filePath,'/instrument_parameters/time_hhmmss').';
        %  GPS time packed,
        % Example: 153320.100 = 15 hours 33 minutes 20 seconds 100 milliseconds.
        
end
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% END OF CODE
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end