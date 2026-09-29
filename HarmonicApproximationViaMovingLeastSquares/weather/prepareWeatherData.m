%% Weather Data: Preparation
% This script writes the two data files used by |weatherExample.m|:
%
%  data/Weather_Data.mat - station positions |nodes| (@vector3d) and daily
%                          mean temperatures |val| in degrees Celsius
%  data/coastLines.mat   - coast lines |coastLines| (@vector3d) for the plots
%
% Both files are included, so this script is only needed to change the date,
% the measured quantity or the coast lines. It downloads the yearly archive
% of the Global Historical Climatology Network - Daily (GHCN-Daily), which is
% about 1 GB, extracts a single day and deletes the archive afterwards. The
% coast lines are the Natural Earth 1:50m land polygons in
% |data/land_contours|; reading them requires the Mapping Toolbox. The script
% also needs |gunzip| and has to be run from the folder |weather|.
%
% Data source: M.J. Menne et al., Global Historical Climatology Network -
% Daily (GHCN-Daily), Version 3, NOAA National Centers for Environmental
% Information, <https://doi.org/10.7289/V5D21VHZ>.

clear
close all

saveVars = true;      % save the results
preview = true;       % plot the extracted data

doDownload = true;    % download the archive and the station list
doExtract = true;     % extract one day -> Weather_Data.mat
doCoastLines = true;  % read the shapefile -> coastLines.mat
doCleanup = true;     % delete the downloaded files

desiredDate = datetime(2018,03,21);  % Fourier's 250th birthday
element = "TAVG";                    % daily mean temperature in tenths of degrees
yr = year(desiredDate);              % the archive is organized by year

dataDir = fullfile(pwd,'data');

url = struct( ...
  'archive',sprintf('https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_year/%i.csv.gz',yr), ...
  'stations','https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt');

file = struct( ...
  'gz',fullfile(dataDir,sprintf('%i.csv.gz',yr)), ...
  'archive',fullfile(dataDir,sprintf('%i.csv',yr)), ...
  'stations',fullfile(dataDir,'ghcnd-stations.txt'), ...
  'shape',fullfile(dataDir,'land_contours','ne_50m_land.shp'), ...
  'data',fullfile(dataDir,'Weather_Data.mat'), ...
  'coastLines',fullfile(dataDir,'coastLines.mat'));

%% Download the archive and the station list

if doDownload
  fprintf('downloading the GHCN archive of %i\n',yr);
  websave(file.gz,url.archive);

  % gunzip replaces the .gz file by the extracted archive
  status = system(sprintf('gunzip -f "%s"',file.gz));
  assert(status == 0,'gunzip failed on %s',file.gz);

  fprintf('downloading the station list\n');
  websave(file.stations,url.stations);
end

%% Extract the measurements of one day
% The archive contains one line per station, day and measured quantity. It is
% sorted by date, so we read it in chunks and stop after the desired date.

if doExtract
  % station positions from the fixed width columns of the station list
  fid = fopen(file.stations,'r');
  C = textscan(fid,'%11s %8f %9f %6f %*[^\n]');
  fclose(fid);
  stationID = string(C{1});
  latS = C{2};
  lonS = C{3};

  % measurements
  ds = tabularTextDatastore(file.archive,'Delimiter',',','ReadVariableNames',false);
  ds.TextscanFormats = {'%s';'%s';'%s';'%f';'%s';'%s';'%s';'%s'};
  ds.VariableNames = {'STATION','DATE','ELEMENT','DATA_VALUE', ...
    'MFLAG','QFLAG','SFLAG','OBSTIME'};
  ds.SelectedVariableNames = {'STATION','DATE','ELEMENT','DATA_VALUE'};

  lon = [];
  lat = [];
  temp = [];

  while hasdata(ds)
    T = read(ds);

    T = T(ismember(T.ELEMENT,element),:);
    if isempty(T), continue, end

    T.DATE = datetime(T.DATE,'InputFormat','yyyyMMdd');
    if T.DATE(1) > desiredDate, break, end

    T = T(T.DATE == desiredDate,:);
    if isempty(T), continue, end

    % tenths of degrees Celsius -> degrees Celsius
    T.DATA_VALUE = T.DATA_VALUE / 10;

    % one row per station, averaging repeated measurements
    P = unstack(T,'DATA_VALUE','ELEMENT','AggregationFunction',@mean);
    if ~all(ismember(element,P.Properties.VariableNames)), continue, end

    % attach the station positions, dropping unlisted stations
    [known,id] = ismember(P.STATION,stationID);
    P = P(known,:);
    id = id(known);

    lon = [lon; lonS(id)];      %#ok<AGROW>
    lat = [lat; latS(id)];      %#ok<AGROW>
    temp = [temp; P.(element)]; %#ok<AGROW>
  end

  fprintf('%i stations\n',numel(temp));

  nodes = vector3d.byPolar(pi/2 - lat*degree,lon*degree);
  nodes.how2plot.outOfScreen = xvector;
  nodes.how2plot.east = yvector;
  val = temp;

  if saveVars, save(file.data,'nodes','val'); end

  if preview
    figure
    plot(nodes,val,'upper','nolabel');
    colormap(WhiteJetColorMap);
    mtexColorbar;
  end
end

%% Coast lines
% The vertices of all polygons are stored in a single @vector3d; the plots
% draw them as small black dots.

if doCoastLines
  shape = shaperead(file.shape,'UseGeoCoords',true);

  rho = cell2mat({shape.Lon})';
  theta = 90 - cell2mat({shape.Lat})';

  % the polygons are separated by NaN
  known = ~(isnan(rho) | isnan(theta));

  coastLines = vector3d.byPolar(theta(known)*degree,rho(known)*degree);

  if saveVars, save(file.coastLines,'coastLines'); end
end

%% Delete the downloaded files

if doCleanup
  for f = {file.archive, file.stations, file.gz}
    if isfile(f{1}), delete(f{1}); end
  end
end
