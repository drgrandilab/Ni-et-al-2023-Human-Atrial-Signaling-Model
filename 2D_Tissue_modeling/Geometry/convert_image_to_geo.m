close all;
d = imread('cell_hetero.png');
figure (1)
imshow(d);





cmap = colormap('jet');
cmap = [0,0,0;cmap];
% cmap = cmap(1:end-8, :);
a = rgb2ind(d, cmap, 'nodither');
% a(a<20)=0;

figure,imshow(a,[])

% J = imresize(a, 0.8);
% imshow(J,[])
% sum(J>0)
% sum(sum(J>0))
% J = imresize(a, 0.7);
% sum(sum(J>0))
% J = imresize(a, 0.6);
% sum(sum(J>0))
% imshow(J,[])
J = imresize(a, 0.62);


% J((J<20) & (J>0) )=20;
J(J<20)=0;
figure, imshow(J,[])
sum(sum(J>0))
geo= J;
% geo (geo>0) = 1;

% dlmwrite('geo.dat', geo(:), 'delimiter', ' ');

cmap = colormap('jet');


cmap = [0,0,0;cmap];
act = imread('activation.png');
figure (4)
imshow(act);
b = rgb2ind(act, cmap, 'nodither');
figure (5)
imshow(b,[]);



x = b;
x(a==0)=0;  % activation map mapped to the geometry (constructed from APD map) map
figure(6);
imshow(x,[])

J = imresize(x, 0.62);
J(J~=8)=0;

J = imresize(x, 0.62);
J(J~=8)=0;
imshow(J,[])
J(8,83)
J(83,8)
J(83,8)=0
imshow(J,[])
J(8,74)=0
imshow(J,[])
J(38,99)=0
imshow(J,[])
J(99,38)=0
imshow(J,[])
J(98,38:42)=0
imshow(J,[])
J(96,38:42)=0
imshow(J,[])
J(J==8)=1;
imshow(J,[])
% dlmwrite('stimulationmap.dat', J(:), 'delimiter', ' ');


[X,Y] = size(J)

newJ = J;
for i = 1:X
    for j = 1:Y
        if( (i-75)*(i-75) + (j-75)*(j-75) > 25)
            newJ(i,j) =0;
        end
    end
end
figure, imshow(newJ,[])
dlmwrite('stimulationmap_small.dat', newJ(:), 'delimiter', ' ');