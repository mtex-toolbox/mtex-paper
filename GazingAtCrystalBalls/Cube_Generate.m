function [screen_int,facedata] = Cube_Generate(BinFile,isHex)
%MasterCube - Read a bin file from Dynamics & Build Interpolants

fileID = fopen(BinFile, 'r', 'ieee-le');
% if fileID == -1, error('Cannot open file: %s', filename); end
format = 'uint';
Data = fread(fileID, Inf, format);
fclose(fileID);
%find out the simulation resolution
cube_res=sqrt((size(Data,1)-12)/6)-1;

%build the face data for the cube
fd=zeros(cube_res+1,cube_res+1,6);
for n=1:6
    dstart=(n-1)*(cube_res+1)^2+1+2*n-1;
    dend=n*(cube_res+1)^2+2*n-1;
    fd(:,:,n)=flipud(rot90(reshape(Data(dstart:dend),cube_res+1,cube_res+1)));
end

%normalise
fd=fd-min(fd(:));
fd=fd/max(fd(:));
facedata=fd;
%sort this data
facedata(:,:,1)=fd(:,:,3); %x+
facedata(:,:,2)=fd(:,:,5); %y+
facedata(:,:,3)=fd(:,:,1); %z+
if isHex == 1 %add in a hexagonal fix
    facedata(:,:,4)=rot90(fd(:,:,1),2);
    facedata(:,:,5)=rot90(fd(:,:,2),2);
    facedata(:,:,6)=(fd(:,:,3));
else
    facedata(:,:,4)=fd(:,:,4); %x-
    facedata(:,:,5)=fd(:,:,6); %y-
    facedata(:,:,6)=fd(:,:,2); %z-
end


%build the interpolants
[gx,gy]=ndgrid(linspace(-1,1,cube_res+1),linspace(-1,1,cube_res+1));
screen_int.p1=griddedInterpolant(gx,gy,facedata(:,:,1),'cubic'); %x+
screen_int.p2=griddedInterpolant(gx,gy,facedata(:,:,2),'cubic'); %y+
screen_int.p3=griddedInterpolant(gx,gy,facedata(:,:,3),'cubic'); %z+
screen_int.p4=griddedInterpolant(gx,gy,facedata(:,:,4),'cubic'); %x-
screen_int.p5=griddedInterpolant(gx,gy,facedata(:,:,5),'cubic'); %y-
screen_int.p6=griddedInterpolant(gx,gy,facedata(:,:,6),'cubic'); %z-
  


end

