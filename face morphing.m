clear;
clc;
close all;
im1 = imread('1.jpg');
im2 = imread('2.jpg');
[H,W,C]=size(im1);
im2=imresize(im2,[H,W]);
%% 
%select correspondence points
cpselect(im1,im2);
%%
load data;
cornerPoints=[1,1;W,1;1,H;W,H];
movingPoints=[movingPoints;cornerPoints];
fixedPoints=[fixedPoints;cornerPoints];
%compute a delauney triangulation based on the correspondences
tri=delaunay(movingPoints);
%%
% aviobj =avifile('morphing.avi');
k=0;
for alpha=0:0.05:1
    k=k+1;
    intermPoints=(1-alpha)*movingPoints+alpha*fixedPoints;

%     hold off;
%     figure(1);imshow(im1);hold all;
%     figure(1);triplot(tri,movingPoints(:,1),movingPoints(:,2),'r');axis image;
% %     figure(1);triplot(tri,intermPoints(:,1),intermPoints(:,2),'g');hold all;
% % 
%     figure(2);imshow(im2);hold all;
%     figure(2);triplot(tri,fixedPoints(:,1),fixedPoints(:,2),'b');hold all;
%     figure(2);triplot(tri,intermPoints(:,1),intermPoints(:,2),'g');

%1 2 3 represent the first one, the second one, the output one respectively
%mapping each triangle
    map1=zeros(H,W);
    map2=zeros(H,W);
    map3=zeros(H,W);
    mask1=im2bw(zeros(H,W));
    mask2=im2bw(zeros(H,W));
    mask3=im2bw(zeros(H,W));
    [X,Y]=meshgrid(1:W,1:H);
    for i=1:size(tri,1)
        vertexX=movingPoints( tri(i,:),1 );
        vertexY=movingPoints( tri(i,:),2 );
        curMap=inpolygon(X,Y,vertexX,vertexY);
        map1=map1+curMap*i.*(1-mask1);
        mask1=(map1~=0);

        vertexX=fixedPoints( tri(i,:),1 );
        vertexY=fixedPoints( tri(i,:),2 );
        curMap=inpolygon(X,Y,vertexX,vertexY);
        map2=map2+curMap*i.*(1-mask2);
        mask2=(map2~=0);

        vertexX=intermPoints( tri(i,:),1 );
        vertexY=intermPoints( tri(i,:),2 );
        curMap=inpolygon(X,Y,vertexX,vertexY);
        map3=map3+curMap*i.*(1-mask3);
        mask3=(map3~=0);
    end

    inputWarp=intermPoints-movingPoints;
    baseWarp=intermPoints-fixedPoints;
    ux1=zeros(H,W);
    uy1=zeros(H,W);
    ux2=zeros(H,W);
    uy2=zeros(H,W);
    for i=1:H
        for j=1:W
            curTri=map1(i,j);
            vertexX=movingPoints( tri(curTri,:),1 );
            vertexY=movingPoints( tri(curTri,:),2 );
            dist=sqrt((vertexX-j).^2+(vertexY-i).^2)+1e-12;
            weight=(1./dist)/sum(1./dist);
            ux1(i,j)=sum(weight.*inputWarp( tri(curTri,:),1 ));
            uy1(i,j)=sum(weight.*inputWarp( tri(curTri,:),2 ));

            curTri=map2(i,j);
            vertexX=fixedPoints( tri(curTri,:),1 );
            vertexY=fixedPoints( tri(curTri,:),2 );
            dist=sqrt((vertexX-j).^2+(vertexY-i).^2)+1e-12;
            weight=(1./dist)/sum(1./dist);
            ux2(i,j)=sum(weight.*baseWarp( tri(curTri,:),1 ));
            uy2(i,j)=sum(weight.*baseWarp( tri(curTri,:),2 ));
        end
    end

    im1Warp=zeros(H,W,C);
    im2Warp=zeros(H,W,C);
    for i=1:C
        im1Warp(:,:,i)=interp2(X,Y,double(im1(:,:,i)),X-ux1,Y-uy1,'bicubic');
        im2Warp(:,:,i)=interp2(X,Y,double(im2(:,:,i)),X-ux2,Y-uy2,'bicubic');
    end
    I=find(isnan(im1Warp));
    im1Warp(I)=zeros(size(I));
    I=find(isnan(im2Warp));
    im2Warp(I)=zeros(size(I));

    
    imBlend=(1-alpha)*im1Warp+alpha*im2Warp;
    imBlend=uint8(imBlend);

    figure(1);imshow(imBlend);
%     aviobj =addframe(aviobj,imBlend);
end
% aviobj =close(aviobj);