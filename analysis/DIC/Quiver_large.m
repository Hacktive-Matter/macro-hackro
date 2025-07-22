clc, clear

% Scale = 124.8305 pixels/mm
% 0.7mm = 87.3813 pixels

A = importdata('RFvalidx.dat'); % x position of raster points in m by n
B = importdata('RFvalidy.dat'); % y position of raster points in m by n

sizevalidx=size(A);
validxfirst=zeros(size(A));
validxfirst=mean(A(:,1),2)*ones(1,sizevalidx(1,2)); % x position of raster points in the reference image

sizevalidy=size(B);
validyfirst=zeros(size(B));
validyfirst=mean(B(:,1),2)*ones(1,sizevalidy(1,2)); % y position of raster points in the reference image

W = dlmread('GridRes.txt');
R = dlmread('RF.txt');
xn = W(1)/R; % reduced resolution in x
yn = W(2)/R; % reduced resolution in y

% we can consider the image with the coordinate system where the upper left corner is its origin.
% displacement of the rasterpoint at the upper left (e.g., @ the coordinate system of (1,1))
xmin = validxfirst(1,1); 
ymin = validyfirst(1,1);

% displacement of the rasterpoint at the bottom right (e.g., @ the coordinate system of (m,n))  
xmax = validxfirst(sizevalidx(1),sizevalidx(2)); 
ymax = validyfirst(sizevalidy(1),sizevalidy(2));

% the dimension of the matrix of the raster points. Use it to reshape obtained displacement which is in m by 1 matrix.
rx = ceil((xmax - (xmin-1))/xn); 
ry = ceil((ymax - (ymin-1))/yn);

% writerObj = VideoWriter('DIC_3fps.mp4');
%  writerObj.FrameRate = 3;
% 
%  % open the video writer
%  open(writerObj);
% 

for i = 1:34
    
    InProgress = i
    second = i+1;
    first = 1;
    
    refx = A(:,first); % reference x position
    comx = A(:,second); % comparing x position
    refy = B(:,first); % reference y position
    comy = B(:,second); % comparing y position
    
    dx = (comx-refx)*(9/124.8305); % x-displacement multiplied by the reduction factor and converted into the length (mm) from pixel 
    dy = (comy-refy)*(9/124.8305); % y-displacement multiplied by the reduction factor and converted into the length (mm) from pixel
    
    x = linspace(1,rx,rx)*W(1)/124.8305; % actual size of the image
    y = linspace(1,ry,ry)*W(1)/124.8305; % actual size of the image
    [X,Y] = meshgrid(x,y);
    u = reshape(dx,ry,rx)*10; % reson for multiplication of 10 is to increase the size of the arrows. Instead, the color bar is reduced by 10 to cancel out.
    v = reshape(dy,ry,rx)*10; % reson for multiplication of 10 is to increase the size of the arrows. Instead, the color bar is reduced by 10 to cancel out.

    mag = sqrt(u.^2+v.^2);
       
    figure(i)
    Q = quiverwcolorbar(X',Y',u',v');
    title(['Displacement Field for Image ',num2str(first),'\rightarrow',num2str(second)])
    set(gca,'FontSize',15)
    xlim([-4 4]); ylim([-1 25])
    
    MAX(:,i) = max(mag,[],"all"); % maximum deformation
    meandef(:,i) = mean(mag,"all"); % average deformation

    % F = getframe(fig);
    % 
    % writeVideo(writerObj,F)
 
end

 % close(writerObj)

% prestr = [6.36 7.06 7.76 8.46 9.16 9.86 10.56 11.26];
% MAX
% meandef
% strain = meandef./prestr

colormin = 0;  
% colormax = 4; % correlation with neighboring image (i.e. relative strain)
colormax = 25; % correlation with reference image (i.e. total strain)

figure(100)
rang = (colormax-colormin)/colormax;
ticknum = 6;        %if you want to toggle number of ticks on colorbar
incr = rang./(ticknum-1);
B = [colormin/colormax:incr:1];
B = B.*colormax;
C = sprintf(['%4.2e',repmat([' \n%4.2e'], 1, ticknum)],B);
C = str2num(C)/10;
caxis([colormin colormax])
colorbar('ytick',B,'yticklabel',C,'ticklength',[0.04 0.1],'YLim',[B(1) B(ticknum)],'FontSize',20)

hold on;
cmap = jet(64);     %toggle type of colormap
CC = colormap(cmap);
cm_stepsize = (colormax-colormin)/length(CC);

% figure(10) % displacement map without color bar
% quiver(X',Y',u',v','LineWidth',1.2,'AutoScaleFactor',2);
% set(gca,'Ydir','reverse','FontSize',13)
% xlim([-3 50]), ylim([-3 38])
% % xlabel('# of grid  in x position'), ylabel('# of grid in y position')

% figure(2), quiverwcolorbar(X',Y',u',v',scale);
% title('Displacement Field for Image 1-2')

% figure(1)
% quiver(X,Y,u1,v1), title(['Total Displacement Field with ', num2str(W(1)),'X',num2str(W(2)), ' grid resolution'])
% xlabel('# of grid  in x position'), ylabel('# of grid in y position')

% figure(1) 
% quiver(X,Y,u2,v2), title(['Average Displacement Field with ', num2str(W(1)),'X',num2str(W(2)), ' grid resolution'])
% xlabel('# of grid  in x position'), ylabel('# of grid in y position')

% % fdx = flip(dx);
% fdx = flip(dx,1);
% % fdy = flip(dy);
% fdy = flip(dy,1);
% m = reshape(fdx,ry,rx);
% n = reshape(fdy,ry,rx);

% figure(2) %flip the graph
% quiver(X,Y,m,n), title(['Flipped Displacement Vector Field with ', num2str(xn),'X',num2str(yn), ' grid resolution'])
% xlabel('# of grid  in x position'), ylabel('# of grid in y position')

% G = importdata('grid_x.dat')

% figure(3), contour(Z,8,'Fill','off','ShowText','on'), set(gca,'Ydir','reverse')

% figure(4), heatmap(Z,'GridVisible','off','Colormap',jet)
% 
% T = 1:1:rx;
% T1 = T(:,2:17);
% T2 = T(:,17:27);
% figure(5), plot(T1,HL1,'k'), hold on, plot(T2,HL2,'k')
% xlabel('X position (pixels)'), ylabel('Displ.Mag. (pixels)')
% 
% O = 1:1:ry;
% O1 = O(:,1:10);
% O2 = O(:,10:22);
% figure(6), plot(O1,VL1,'k'), hold on, plot(O2,VL2,'k')
% xlabel('Y position (pixels)'), ylabel('Displ.Mag. (pixels)')













