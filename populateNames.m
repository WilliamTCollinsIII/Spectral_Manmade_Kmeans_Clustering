% Custom update function for data cursor
function txt = populateNames(~, event_obj, manmadeNames, centroids, centroidClusterAssignments)
    % Gets the index of the selected data point
    idx = get(event_obj, 'DataIndex');
    
    % Retrieves the position of the selected data point
    pos = get(event_obj, 'Position');
    
    % Check if the selected point is a centroid
    [isCentroid, centroidIdx] = ismember(pos, centroids, 'rows');
    
    if isCentroid
        % If the point is a centroid, display the cluster number
        txt = {['Cluster: ', num2str(centroidClusterAssignments(centroidIdx))]};
    else
        % If not a centroid, return the name
        name = manmadeNames{idx};
        txt = {['Name: ', name]};
    end
end