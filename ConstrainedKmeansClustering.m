%{

Title: Exact Name Constrained Kmeans Clustering
Author: William Collins

The purpose of this script is to:
- Read in the ASTER and KLUM data sets.
- Down sample the data to the wavelengths desired
- Group samples together that have the same file name
- Perform PCA to reduce the dimensionality of the data
- Perform constrained Kmeans clustering so the groups are always assigned
the same cluster.

%}

clc
clear
close all

path = pwd;
temp = dir([path, '\*.spectrum.txt']);

% Create cell arrays for each file group
manmadeFiles = {};

% Set counts to 0 for each file group
manmadeCount = 0;

% Read in file names
fileNames = {temp.name};

% Read files into respective cell array groups
for i = 1:length(fileNames)
    if contains(fileNames{i}, 'manmade')
        manmadeFiles{end+1} = fileNames{i};
        manmadeCount = manmadeCount + 1;
    end
end

% Read KLUM dataset
[klumWavelength, klumReflectivity, klumNames] = readKlum();
klumCount = size(klumReflectivity,1);

% Create and pre-allocate cell arrays for wavelength and reflectivity data
manmadeWavelength = cell(1, manmadeCount);
manmadeReflectivity = cell(1, manmadeCount);

% Reading and Filtering the data
[manmadeBadCount, manmadeCount, manmadeWavelength, manmadeReflectivity, asterNames] = filterAndRead(manmadeCount, manmadeFiles, manmadeWavelength, manmadeReflectivity);
asterNames = asterNames';

manmadeNames = [asterNames; klumNames];

% Add KLUM to manmade
for i = manmadeCount + 1:manmadeCount + klumCount
    manmadeWavelength{i} = klumWavelength(i-manmadeCount,:)';
    manmadeReflectivity{i} = klumReflectivity(i-manmadeCount,:)';
end

identity = 1;
[manmadeReflectivityMatrix, manmadeWavelengthMatrix, manmadeIdentityMatrix] = knnMatrixSetUp(manmadeCount + klumCount, manmadeReflectivity, manmadeWavelength, identity);

% Downsample to wavelengths (change values if desired)
wavelengths = {1.06 1.08 1.10 1.22 1.24 1.26 1.28 1.30 1.52 1.54 1.56 1.58 1.60 1.62 1.64}; 
[manmadeWavelengthMatrix, manmadeReflectivityMatrix] = downsample2Wavelengths(manmadeWavelengthMatrix, wavelengths, manmadeReflectivityMatrix);

% Normalize the data
manmadeReflectivityMatrix = cell2mat(manmadeReflectivityMatrix);
normalizedManmade = normalize(manmadeReflectivityMatrix);
manmadeReflectivityMatrix = num2cell(normalizedManmade);

% Put all data in one matrix
dataMatrix = manmadeReflectivityMatrix;
identityMatrix = manmadeIdentityMatrix;
wavelengthMatrix = manmadeWavelengthMatrix;

[dataMatrix, wavelengthMatrix, identityMatrix] = removeNans(dataMatrix, wavelengthMatrix, identityMatrix);
[wavelengthMatrix, dataMatrix] = sortData(dataMatrix, wavelengthMatrix);

X = cell2mat(dataMatrix);

X_standardized = zscore(X);

numClusters = 10;
%for numClusters = 8: 27
% Initialize an empty table for cluster information
clusterInfo = table();

% Specify the path for the Excel file
excelFilePath = fullfile(path, 'myClusterAssignments.xlsx');

% Perform PCA with numComponents
[coeff, score, ~, ~, explained] = pca(X_standardized, 'NumComponents', 3);
reducedData = score(:, 1:3);

% Implement constrained K-means clustering
maxIter = 500; % Maximum number of iterations
tol = 1e-6; % Tolerance for convergence
numPoints = size(reducedData, 1);

% Randomly initialize the centroids
rng(1); % For reproducibility
centroids = reducedData(randperm(numPoints, numClusters), :);

% Initialize variables
clusterAssignments = zeros(numPoints, 1);
previousCentroids = centroids;

% Group indices by names
uniqueNames = unique(manmadeNames);
nameGroups = cell(length(uniqueNames), 1);
for i = 1:length(uniqueNames)
    nameGroups{i} = find(strcmp(manmadeNames, uniqueNames{i}));
end

for iter = 1:maxIter

   for groupIdx = 1:length(nameGroups)
    group = nameGroups{groupIdx};
    % Calculate distances for each point in the group to all centroids
    distances = zeros(length(group), numClusters);
        for j = 1:numClusters
            distances(:, j) = sum((reducedData(group, :) - centroids(j, :)).^2, 2);
        end
    % Assign the entire group to the nearest centroid
    [~, minIndex] = min(distances, [], 2);
    clusterAssignments(group) = mode(minIndex); % Assign the mode of minIndex to ensure group consistency

    end
    % Update step: Calculate new centroids
    previousCentroids = centroids;
    for j = 1:numClusters
        pointsInCluster = reducedData(clusterAssignments == j, :);
        if ~isempty(pointsInCluster)
            centroids(j, :) = mean(pointsInCluster, 1);
        end
    end

    % Check for convergence
    if max(max(abs(centroids - previousCentroids))) < tol
        break;
    end
end

% Create a table with file names and cluster assignments
clusterAssignmentsTable = table(manmadeNames, clusterAssignments, 'VariableNames', {'FileName', 'ClusterAssignment'});

% Sort the table by cluster assignments
sortedClusterAssignments = sortrows(clusterAssignmentsTable, 'ClusterAssignment');

% Generate a sheet name for this specific cluster and component configuration
sheetName = sprintf('%dAssignments',numClusters);

% Write the table to the specified sheet in the Excel file
writetable(sortedClusterAssignments, excelFilePath, 'Sheet', sheetName);
%end

% Create variable to populate centroid number
centroidClusterAssignments = (1:numClusters)';

% Plot the clustered data
figure;
scatter3(reducedData(:, 1), reducedData(:, 2), reducedData(:, 3), 36, clusterAssignments, 'filled');
% Enable data cursor mode
dcm_obj = datacursormode(gcf);
set(dcm_obj, 'UpdateFcn', {@populateNames, manmadeNames, centroids, centroidClusterAssignments});
hold on;
scatter3(centroids(:, 1), centroids(:, 2), centroids(:, 3), 100, 'k', 'filled', 'MarkerEdgeColor', 'w');
title('Constrained K-means Clustering on PCA-reduced Data');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
zlabel('Principal Component 3');
colormap jet;
hold off;