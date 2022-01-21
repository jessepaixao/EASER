% Running SLICER
% File I/O
tic
clear triangles movelist
triangles = read_stl_file(filename);
timeio=toc;
fprintf('File read done, %.4f sec elapsed\n',timeio);
% Fratio,
% Rotate Model
triangles = rotate_model(triangles,'x',0);
triangles = [triangles(:,1:12),min(triangles(:,[3 6 9]),[],2), max(triangles(:,[ 3 6 9]),[],2)];
timemr=toc;
fprintf('Model rotate done, %.4f sec elapsed\n',timemr);

% Generate layers
if adptive==1
    [movelist(:,1), movelist(:,2), movelist(:,3)] = adptive_polygon(triangles);
else
    [movelist(:,1), movelist(:,2), movelist(:,3)] = generate_polygon(triangles);
end
timegl=toc;
fprintf('Generate layer done, %.4f sec elapsed\n',timegl);

% Generate shell
[movelist(:,4)] = generate_shell(movelist(:,3));
timegs=toc;
fprintf('Generate shell done, %.4f sec elapsed\n',timegs);

% Generate infill
x_min=min(min(triangles(:,1:3:9)));
x_max=max(max(triangles(:,1:3:9)));
y_min=min(min(triangles(:,2:3:9)));
y_max=max(max(triangles(:,2:3:9)));
row_index=find([movelist{:,1}]==1);
[movelist(row_index,5)] = generate_infill(movelist(row_index,3), x_min, x_max, y_min, y_max);
% close all;
timegi=toc;
fprintf('Generate infill done, %.4f sec elapsed\n',timegi);

% Export gcode
statu=exp_gcode2(movelist(:,2), movelist(:,4), movelist(:,5), gfilename);
timeeg=toc;
if statu==1
    fprintf('Gcode file has written to %s\nElapsed time %.4f sec\n',gfilename,timeeg);
end
