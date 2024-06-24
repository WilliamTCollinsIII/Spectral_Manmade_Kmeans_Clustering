function [dataNames, newReflectivityData, newWavelengthData] = readArmyData()
clc;
clear;

% Specify the file path
filePath = "C:\Users\willc\OneDrive\Documents\AFRL\Kmeans\NVL_ARMY_Paints.xlsx";

% Load the "paint" sheet
opts = detectImportOptions(filePath, 'Sheet', 'paint');

% Read the data from the sheet
dataTable = readtable(filePath, opts, 'Sheet', 'paint');

% Find all the columns that contain "Wavelength"
wavelengthColumns = contains(dataTable.Properties.VariableNames, 'Wavelength');

% Extract the wavelength values (assuming they are in every third column starting from 1)
wavelengthsData = dataTable{:, wavelengthColumns};

% Extract the name of data (column headers excluding "Wavelength")
dataNames = dataTable.Properties.VariableNames(~wavelengthColumns)';

% Extract the reflectivity data (all columns excluding the ones with "Wavelength")
reflectivityData = dataTable{:, ~wavelengthColumns};

% Specify the new wavelengths for which we need to interpolate the data
Wavelengths = [1.06 1.08 1.10 1.22 1.24 1.26 1.28 1.30 1.52 1.54 1.56 1.58 1.60 1.62 1.64];

% Initialize the array to store the interpolated reflectivity data
newReflectivityData = zeros(size(reflectivityData, 2), length(Wavelengths));

% Loop over each reflectivity dataset and interpolate
for i = 1:size(reflectivityData, 2)

    % Get the current wavelength and reflectivity data
    currentWavelengths = wavelengthsData(:, i);
    currentReflectivity = reflectivityData(:, i);
    
    % Perform the interpolation
    newReflectivityData(i, :) = interp1(currentWavelengths, currentReflectivity, Wavelengths, 'linear');
    newWavelengthData(i,:) = Wavelengths;
    for j = 1 : length(Wavelengths)
        newReflectivityData(i,j) = newReflectivityData(i,j)*100;
    end
end
end