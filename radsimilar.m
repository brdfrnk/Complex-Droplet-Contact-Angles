% Contact Angle Attempt One
%   This version finds all curves, does some math to attach two closest
%   curves to attempt to find droplet contact angle.

%   On opening an image you will be asked to draw a box over the radius of
%   one droplet.

% Data gathering and image loading
clc; clear all; close all;

r_min = 0.8;
r_max = 1.8;
sensitivity = 0.97;

[filename,path] = uigetfile({'*.jpg';'*.png'},'Select Image for Analysis');
imageinput=strcat(path, filename);
image = imread(imageinput);
imshow(image)

hold on
title('Please trace the radii of one droplet')
[height,width,colors] = size(image);
set(gcf,'Units','normalized')
k = waitforbuttonpress; %Ask user to draw a box
rect = rbbox;
rad = (rect(1,3)*width)/2;

radius_min = round(((rad))*r_min);
radius_max = round(((rad))*r_max);
imshow(image);
hold on
title(filename)
[centers,radii] = imfindcircles(image,[radius_min radius_max], 'ObjectPolarity','dark','Sensitivity',sensitivity);

circle_neighbors = rangesearch(centers,centers,radius_max); %searches each center for a neighbor, using half the minimum radius
for i = 1:size(circle_neighbors) %Sorting out only pairs, not singles, not triplets
    if length(circle_neighbors{i}) == 2
        %its good
    else
        circle_neighbors{i} = [];
    end
end
circle_pairs = circle_neighbors(~cellfun('isempty',circle_neighbors)); %Remove everything with >2 entries
for i = 1:size(circle_pairs) %Sort the entries by row in order to remove duplicate pairs
    circle_pairs{i} = sort(circle_pairs{i},'ascend');
end
unique_pairs = unique(cell2mat(circle_pairs),'rows'); %Finally, only keep the unique pairs each in a row

if length(unique_pairs) >= 1
    unique_total = length(unique_pairs(:,1));
    pairs_a = zeros(unique_total,2);
    pairs_b = zeros(unique_total,2);
    radii_a = zeros(unique_total,1);
    radii_b = zeros(unique_total,1);
    
    for i = 1:unique_total %Get only our 'successful' pairs
        pairs_a(i,:) = [centers(unique_pairs(i,1),1);centers(unique_pairs(i,1),2)];
        pairs_b(i,:) = [centers(unique_pairs(i,2),1);centers(unique_pairs(i,2),2)];
        radii_a(i,:) = radii(unique_pairs(i,1),:);
        radii_b(i,:) = radii(unique_pairs(i,2),:);
    end
    if ~isempty(pairs_a)
        viscircles(pairs_a,radii_a,'Color','w','LineWidth',0.4);
        viscircles(pairs_b,radii_b,'Color','cyan','LineWidth',0.4);
    end
    
    for i=1:unique_total %Displays the # of all detected pairs
        text(pairs_a(i,1),pairs_a(i,2),join([num2str(i),"a"]),'Color','white');
        text(pairs_b(i,1),pairs_b(i,2),join([num2str(i),"b"]),'Color','cyan');
    end
    
    %Actual math section
    angles = zeros(unique_total,3);
    
    for i = 1:unique_total
        midpoint_a = [centers(unique_pairs(i,1),1);centers(unique_pairs(i,1),2)];
        midpoint_b = [centers(unique_pairs(i,2),1);centers(unique_pairs(i,2),2)];
        radius_a = radii(unique_pairs(i,1),:);
        radius_b = radii(unique_pairs(i,2),:);
        
        d2 = sum((midpoint_b-midpoint_a).^2);
        P0 = (midpoint_a+midpoint_b)/2+(radius_a^2-radius_b^2)/d2/2*(midpoint_b-midpoint_a);
        t = ((radius_a+radius_b)^2-d2)*(d2-(radius_b-radius_a)^2);
        if t <= 0
            %fprintf('No intersection\n')
        else
            T = sqrt(t)/d2/2*[0 -1;1 0]*(midpoint_b-midpoint_a);
            intersect = P0 - T;
        end
        if t > 0
            angle = (atan2(abs((midpoint_a(1)-intersect(1))*(midpoint_b(2)-intersect(2))-(midpoint_b(1)-intersect(1))*(midpoint_a(2)-intersect(2))), ...
                (midpoint_a(1)-intersect(1))*(midpoint_b(1)-intersect(1))+(midpoint_a(2)-intersect(2))*(midpoint_b(2)-intersect(2))))*180/3.14159;
            both_angle_cases = [angle 180-angle];
            
            % text(midpoint_a(1,1),midpoint_a(2,1),num2str(unique_pairs(i,1)),'Color','white');
            % text(midpoint_b(1,1),midpoint_b(2,1),num2str(unique_pairs(i,2)),'Color','cyan');
            scatter(intersect(1,1),intersect(2,1))
            text(intersect(1,1),intersect(2,1),num2str(both_angle_cases(1,1)),'Color','red');
            angles(i,:) = [i, angle, 180-angle];
        end
    end
    hold off
    csvout = strcat('output',filename(1:end-4), '-out.csv');
    pngout = strcat('output',filename(1:end-4), '-out.png');
    writematrix(angles,csvout);
    export = gcf;
    saveas(export,pngout);
else
    unique_total = 0;
end
total_pairs = unique_total