function [statu]= exp_gcode2(height,shell,infill,gfilename)
global bed_temp filament_temp nozzle_dim layer_hight filament_dim print_speed init_z bed_x0 bed_z0
% Eratio=nozzle_dim*layer_hight/pi/(filament_dim/2)^2;
Fratio=print_speed*30;
statu=0;
fid = fopen(gfilename,'w');
% fprintf(fid,'M104 S%d ;extruder temperature\n',filament_temp);
% fprintf(fid,'M109 S%d ;extruder temperature, and wait\n',filament_temp);
fprintf(fid,'M82 ;absolute extrusion\n');
fprintf(fid,'G28 ;home\n');
fprintf(fid,'G1 Z%d F3000.000\n',init_z);
fprintf(fid,'G92 E10; \n');
% fprintf(fid,'G1 F200 E2\n');
% fprintf(fid,'G92 E0\n');
fprintf(fid,'M83 ;relative extrusion\n');
fprintf(fid,'M107 ;fan off\n');
fprintf(fid,';start\n');
% fprintf(fid,'M106 S255 ;fan on\n');
fprintf(fid,'M104 S%d ;extruder temperature\n',filament_temp);
fprintf(fid,'M109 S%d ;extruder temperature, and wait\n',filament_temp);
% fprintf(fid,'M104 S%d ;extruder temperature\n',filament_temp);
% fprintf(fid,'G1 E10\n');
% fprintf(fid,'M140 S%d ;bed temperature\n',bed_temp);
gstatu=2;%0 normal,1 try to fix,2 end of loop or end of layer,3 end of shell
for i=1:size(height,1)
    Z=height{i};
    point=[shell{i};NaN,NaN;infill{i}];
    if i==1
        Eratio=3*nozzle_dim*layer_hight/pi/(filament_dim/2)^2;
%         Eratio=0.5;
    else
        Eratio=3.5*nozzle_dim*(height{i}-height{i-1})/pi/(filament_dim/2)^2;
    end
    
    if i==size(height,1)
        Eratio=3*nozzle_dim*layer_hight/pi/(filament_dim/2)^2;
    end
    
    for j=1:size(point,1)
        point(j,1);
        point(j,2);
        if isnan(point(j,1))
            %nan denote end of loop or end of layer
            gstatu=2;
            continue;%continue end current for loop,break end whole for loop
        end
        if  gstatu==1
            x1=point(j,1);y1=point(j,2);
            x0=point(j-1,1);y0=point(j-1,2);
            distance=sqrt((x1-x0)^2+(y1-y0)^2); 
            fprintf(fid,'G1 X%.4f Y%.4f E%.5f\n',point(j,1),point(j,2),distance*Eratio);
            gstatu=0;
            continue;
        elseif gstatu==2 
            fprintf(fid,'G1 X%.4f Y%.4f Z%.4f F%d\n',point(j,1),point(j,2),Z,Fratio*2);
            gstatu=1;
            continue;
        end
        x1=point(j,1);y1=point(j,2);
        x0=point(j-1,1);y0=point(j-1,2);
        distance=sqrt((x1-x0)^2+(y1-y0)^2); 
        fprintf(fid,'G1 X%.4f Y%.4f E%.5f\n',point(j,1),point(j,2),distance*Eratio);
    end
    fprintf(fid,';LAYER: %d\n',i);
    if i==2
%         fprintf(fid,'M106 S255 ;fan on\n');
%         fprintf(fid,'M104 S%d ;extruder temperature\n',filament_temp);
%         fprintf(fid,'M140 S%d ;bed temperature\n',bed_temp);
    end
    gstatu=2;
end
fprintf(fid,';end\n');
% fprintf(fid,'G1 X100.0000 Z150.000 E0.00; lift nozzle\n');

fprintf(fid,'G1 X%.4f ; lift nozzle\n',bed_x0-20);
fprintf(fid,'G1 Z%.4f ; lift nozzle\n',bed_z0+20);
fprintf(fid,'G1 Y0.0000 F1800; lift nozzle\n');
% fprintf(fid,'G1 F1500 E-0.2 ;retract\n');
% fprintf(fid,'M82 ;absolute extrusion\n');
fprintf(fid,'M107 ;fan off\n');
% fprintf(fid,'M140 S0 ;bed temperature 0\n');
fprintf(fid,'M104 S250 ;extruder temperature 0\n');
% fprintf(fid,'G92 E0 ;reset extrusion distance\n');
fprintf(fid,'G28 X0 Y0 ;home\n');
fprintf(fid,'M84 ;motors off\n');
fprintf(fid,'quit ;\n');
fclose(fid);
fclose all;
%movefile test.txt gfilename;
statu=1;
end