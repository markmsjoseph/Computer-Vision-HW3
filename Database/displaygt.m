% set class = 'S' or 'M' to view Single- or Multiple-symmetry database,
%     respectively
% press the space bar to go to the next image

class = 'S'; % 'S' or 'M'

if class == 'S'
    nimages = 176;
elseif class == 'M'
    nimages = 63;
end

for index = 1:nimages
    im_path = [class sprintf('/I%03d.png',index)];
    gt_path = [class sprintf('/I%03d.mat',index)];

    I = imread(im_path);
    load(gt_path)
    
    for i = 1:length(segments)
        points = segments{i};
        p = points(1,:);
        q = points(2,:);
        npq = norm(q-p);
        v = (q-p)/npq;
        for j = 1:npq
            r = round(p+j*v);
            r(2) = min(max(r(2),2),size(I,1)-1);
            r(1) = min(max(r(1),2),size(I,2)-1);
            I(r(2)-1:r(2)+1,r(1)-1:r(1)+1,:) = 0;
        end
        for j = 1:npq
            r = round(p+j*v);
            I(r(2),r(1),:) = 255;
        end
    end
    
    imshow(I)
    title(sprintf('%d',index))
    pause
end
close all