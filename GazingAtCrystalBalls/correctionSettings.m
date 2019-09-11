%background correctionData_InputMapData_InputMap
Settings_Cor.gfilt=1; %use a low pass filter (do you mean high pass?)
Settings_Cor.gfilt_s=4; %low pass filter sigma

%radius mask
Settings_Cor.radius=1; %use a radius mask
Settings_Cor.radius_frac=0.95; %fraction of the pattern width to use as the mask

%hold pixelData_InputMap
Settings_Cor.hotpixel=0; %hot pixel correction
Settings_Cor.hot_thresh=1000; %hot pixel threshold

%resize
Settings_Cor.resize=1; %resize correction
Settings_Cor.size=150; %image width

Settings_Cor.RealBG=0; %use a real BG
Settings_Cor.EBSP_bgnum=30; %number of real pattern to use for BG

Settings_Cor2=Settings_Cor;
Settings_Cor2.radius=1;
Settings_Cor2.size=300;
Settings_Cor2.SplitBG=1;
Settings_Cor2.gfilt=0;
Settings_Cor2.hotpixel=1;
Settings_Cor2.hot_thresh=1000;