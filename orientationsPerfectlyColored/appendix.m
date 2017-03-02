% Script file for the Appendix of the publication
%
%% IPF coloring of crystal orientation data
%
% Gert Nolze and Ralf Hielscher 
%


%% Appendix A.1

h = linspace(0,1);
s = ones(size(h));
v = ones(size(h));
rgb = reshape(hsv2rgb(h,s,v),1,[],3);
imagesc(rgb), axis off

%saveFigure('hsv_default.png')
%saveFigure('hsv_modified.png')

%% Appendix A.2

cs = crystalSymmetry('1');
oM = ipdfHSVOrientationMapping(cs);
oM.colorPostRotation = rotation('Euler',-pi/2,pi/2,0,'ZYZ');
oM.grayGradient = 1;
oM.grayValue = 0;

plot(oM,'noLabel');
set(gcf,'position',[200 200 230 450])
% saveFigure('../pic/1default.png')

%%

oM.grayGradient = [0.25 0.5];
oM.grayValue = 0;

plot(oM,'noLabel');
set(gcf,'position',[200 200 230 450])

% saveFigure('../pic/1_02.png')

%%

oM.grayGradient = [0.25 0.5];
oM.grayValue = [0.2,0.5];

plot(oM,'noLabel');
set(gcf,'position',[200 200 230 450])

% saveFigure('../pic/1_05_02.png')


%% Appendix B
%
% plot all ipf keys

for id = 1:45

  cs = crystalSymmetry('pointid',id);
  disp(cs.pointGroup)
  oM = ipdfHSVOrientationMapping(cs);

  plot(oM,'noTitle','noLabel','innerplotSpacing',5)
  hold on
  plot(cs)
  %hold on
  %plot(cs.fundamentalSector,'color','red')
  pause
  hold off
  %saveFigure(['../full/' strrep(cs.pointGroup,'/','_') '.png'])
  
end

% 2/m11(ok) 121(ok) 112(ok) -3 312 -31m(ok)



