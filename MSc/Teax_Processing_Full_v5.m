%%%% INPUT %%%
LocIm = 'C:\Users\womaes\Downloads\20240716_Seluleko\20240716_Seluleko\20240716_Thermal\TIF\'; %location of the images

%v2: all images are exported from Thermoviewer, the sharpest per second
%is/are selected; the other images are removed (to save space)
%(Opt_RemUnsh)

%v3: v2 + also working when no images with GPS coordinates are present in the dataset. 

%v4: improved multispectral functionality, scanning in each subfolder + new
%import

%v5% Options: 
% CHANGED AN IMPORTANT ERROR IN THE Ts-FORMULA. ONLY USE V FROM NOW ONE!!
% - Option to output as brightness temperature instead of ts (for
% later NDVI-based emissivity calculation)

Opt_Tcorr = 1;  %put 1 if you want to perform this option. If weather conditions during the flight were highly changing, 
                % use this to correct for Temperature fluctuations which
                % are stored in an additional folder 
Opt_Tbr = 0;   % if 1, you will calculate Brightness temperature, NOT surface temperature. 
                %  This is if you want to use the NDVI-based emissivity correction for more precise Ts estimates.                
Opt_RemUnsh=1; %Unsharp images are removed in the first step of the process 
Opt_Sharp =1;  %If you want to move unsharp images to separate folder
Opt_LocOnly=0; %1 if you want to use the script only to save the location of the cameras; otherwise, 1
Opt_PosMulti=0; %if you want to estimate the XYZ and PitchYawRoll for the thermal cameras based on the multispectral data. 
                % This requires that the multispectral data are processed
                % in METASHAPE, and that you have exported the estimated
                % camera locations. Matching will be done based on the GNSS
                % time. This works only if you have selected the option to
                % create a XXX_meta.csv file when exporting the cameras in
                % Thermoviewer! Regardless, another file is also created
                % with the camera positions recorded by the GPS.

Opt_Csv = 0;    % 1 for csv-data, 0 for Tif data
               
                
% Parameters:

Par.RH = 0;         %Relative Humidity, in percentage. If unknown, try looking it up at https://wow.meteo.be/nl/. Else, put 0
Par.Tair = 0;       %Air temperature, in °C; If unknown, try looking it up at https://wow.meteo.be/nl/. Else, put 0
Par.Tbg = -7.6;        %Background temperature. Seen in Thermoviewer as Tbr of alu panels. 
Par.FH = 17;        %Flight Height - actually height to object (flight height minus crop/tree height)
Par.ImPerSec = 2;   %Images per second you want to keep.

% Only needed for Opt_PosMulti, ie if you want to estimate the XYZ and
% YawPitchRoll of the thermal based on the multispectral data. This works
% only for the TIF option!
LocMulti = 'C:\Users\womaes\Downloads\20240716_Seluleko\20240716_Seluleko\20240716_Multispec\';       % Location of multispectral Tif data. Careful, this needs to be the original data! (Not the reflectance product)
LocTxtMulti = 'C:\Users\womaes\Downloads\20240716_Seluleko\20240716_Seluleko\20240716_Multispec\'; %Location of file with multispectral locationHere
NameTxt = '06_07_24_camera_position_multispec.txt';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra INPUT: don't change unless absolutely necessary and you know what you're doing!
Par.TLimLow = 253; Par.TLimHigh = 353; %in K; Don't change unless absolutely necessary
Par.emiss= 0.985;    %Emissivity
Par.SmoothLength = 60; %Total smoothing for Opt_Tcorr; this is in number of images, with usually 3 images per second. (60= smoothing of 10s before and 10s after point) 
Par.PTileThresh = 3; % Parameters for removing the unsharp data
Par.MaxGap = 4;
Par.StepThresh = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This script reads the data from Thermoviewer (either csv or tif, to be mentioned), performs the full radiometric correction, 
%and saves the Ts-data as Tif files. 

%It then looks at sharp vs unsharp images, and moves less sharp images to a
%separate folder. 

%In case Ts-data were obtained in cloudy/changing conditions, a correction
%(Opt_Tcorr) is added, which performs a smoothing of the Ts-data and
%corrects for temperature fluctuations. This has the danger of also
%eliminating areas with hotter/cooler conditions. A smoothing of 60 frames
%(roughly 20 seconds of flight, or about 60-80m) is now implemented. This
%can be adjusted 

%%%%%%%%%%%%%%% Actual Script
warning ('off','all');

%Create LocExp
if ~(strcmp(LocIm(end),'\') | strcmp(LocIm(end),'/'))
    LocIm = [LocIm '\'];
end
if ~(strcmp(LocTxtMulti(end),'\') | strcmp(LocTxtMulti(end),'/'))
    LocTxtMulti = [LocTxtMulti '\'];
end
if ~(strcmp(LocMulti(end),'\') | strcmp(LocMulti(end),'/'))
    LocMulti = [LocMulti '\'];
end
LocExp = [LocIm 'Tif_Ts\']; mkdir(LocExp);

if Opt_Csv
    DirIm = dir([LocIm '*.csv']);
else
    DirIm = dir([LocIm '*.tif']);
end

%% Import some basic data (time, location) per image; define if NoGPS
DirMeta = dir([LocIm '*_meta.csv']);
if round(length(DirIm)/10,0) == round(length(DirMeta)/10,0) %If by accident you saved a metadata file for each image, this part corrects for this.
   for k1=1:length(DirMeta)
       if k1==1 
           ThermalInf = ImportMeta([LocIm DirMeta(1,1).name]);
       else
           ThermalInf = [ThermalInf; ImportMeta([LocIm DirMeta(k1,1).name])];
       end
   end
else
    ThermalInf = ImportMeta([LocIm DirMeta(1,1).name]);
end
  
clear k1 

 % Check: if no GPS recordings
if ~(any(isnumeric(ThermalInf.Latitude)))
    display('No GPS data recorded');
    NoGPS=1;
else  
    NoGPS=0;
end


%% First loop: only keep the sharpest X images per second (Par.ImPerSec), remove all other
%images

TTime = nan(length(DirIm),1); %Define empty tables first
TSharp = TTime;TTime_corr = nan(size(TTime));KeepT = zeros(size(TTime));

if ~Opt_LocOnly
    for i1=1:length(DirIm)
        A= f_readIm([LocIm DirIm(i1,1).name], Opt_Csv);
        A = double(A);
        [TSharp(i1)]=estimate_sharpness(A);
    end
end
clear A


%% First, get the time of the thermal data capture
if length(DirIm) ~= length(ThermalInf.Date) 
    
    display('Nr of Tifs not the same as in Meta-table');
else
    if ~NoGPS
                    % write down thermal time
                     %Thermal Time
                     DateS=cellstr(ThermalInf.Date(2)); Date = datenum(str2num(DateS{1}(7:10)),str2num(DateS{1}(4:5)),str2num(DateS{1}(1:2)));
                     TTime = Date + datenum(ThermalInf.Time) - floor(datenum(ThermalInf.Time));                 


         %Scan the nr of images taken per sec.
        UnTT = unique(TTime);
         for i7=1:length(UnTT)
           Id = find(UnTT(i7)==TTime);
           for i8 = 1:length(Id)
              TTime_corr(Id(i8)) = TTime(Id(i8)) + (i8-1)/(length(Id).*(24.*3600)); %Exact time for the images, below second level
           end
           [~,KeepId] = maxk(TSharp(Id),Par.ImPerSec);
           KeepT(Id(KeepId))=1;
         end

         clear i1 i7 i8 Id KeepId UnTT Date DateS 

         % Now remove the unsharp images and continue working with the sharp images only.

         TTime_corr = TTime_corr(KeepT==1);
         TTime = TTime(KeepT==1);
         TSharp=TSharp(KeepT==1);
         ThermalInf = ThermalInf(KeepT==1,:);
         RemCam = DirIm(KeepT==0);
         DirIm = DirIm(KeepT==1);

         if Opt_RemUnsh
             for i2=1:length(RemCam)
                 delete([LocIm RemCam(i2).name]);
             end
         end
     
    else
         %If no GPS (or time) stamp, we will group the images per five, and keep the Par.ImPerSec nr of images per group of 5.

           if Opt_RemUnsh
               LDir = round(length(DirIm)/5);
               for i2=1:LDir
                    IDs = [(i2-1)*5+1:1:i2*5];
                    [~,KeepId] = maxk(TSharp(IDs),Par.ImPerSec);
                    for i3=1:Par.ImPerSec
                        KeepT((i2-1)*5+KeepId(i3))=1;
                    end
               end
           end
      
        % Now remove the unsharp images and continue working with the sharp images only.
       
         TSharp=TSharp(KeepT==1);
         ThermalInf = ThermalInf(KeepT==1,:);
         RemCam = DirIm(KeepT==0);
         DirIm = DirIm(KeepT==1);
         
         if Opt_RemUnsh
             for i2=1:length(RemCam)
             %    delete([LocIm Rem(i2).name]);
             end
         end


    end

end
%% 
 clear KeepT RemCam
 %Convert Tbg and ImCsv to K if needed
  if Par.Tbg <200
     Par.Tbg = 273.15+Par.Tbg;
  end

  if Opt_Tcorr
   TAv = nan(size(DirIm));
  end

% calculating transmission
% Code % Obtained from Heinemann et al 2020 RS https://doi.org/10.3390/rs12071075
%Can be placed in the loop, if the Par.FH is made dynamic. 

     if Par.Tair ==0
         A= f_readIm([LocIm DirIm(1,1).name], Opt_Csv);
         Par.Tair = mean(A(:)); 
         if Par.Tair > 200 %convert to K if needed (but shouldn't be!)
             Par.Tair = Par.Tair-273.15;
         end
     end
     
     if Par.RH ==0
         Par.RH=50;
     end

Par.omega = Par.RH/100*exp(6.8455*10^-7*Par.Tair.^3 -2.7816*10^-4*Par.Tair.^2+6.939*10^-2*Par.Tair+1.5587);
Par.Tau_atm = 1.9*exp(-sqrt(Par.FH)*(0.0066 -0.0023*sqrt(Par.omega))) + (1-1.9)*exp(-sqrt(Par.FH)*(0.0126 -0.0067*sqrt(Par.omega)));

 clear A
 
%%
%Loop through image directory: reading thermal images, converting to Ts and
%exporting

if ~Opt_LocOnly

    for i1=1:length(DirIm)
    %  try
        A= f_readIm([LocIm DirIm(i1,1).name], Opt_Csv);
        %% Convert temperature data to K
         if A(1,1)<200
            A = 273.15+A;
         end

         %% Correcting for atmospheric conditions
         if Opt_Tbr
             Ts= ((A.^4 - (1-Par.Tau_atm).*(Par.Tair+273.15).^4)/(Par.Tau_atm)).^0.25;
         else
           
           Ts= ((A.^4 - Par.Tau_atm.*(1-Par.emiss).*Par.Tbg.^4 - (1-Par.Tau_atm).*(Par.Tair+273.15).^4)/(Par.emiss.*Par.Tau_atm)).^0.25;

                   %Convert images to 16-bit tif
          
         end
          TsConv = (Ts -Par.TLimLow)./(Par.TLimHigh-Par.TLimLow).*65535;
           Image16 = uint16(round(TsConv));
           NameIm = [LocExp DirIm(i1,1).name(1:end-4) '.tif'];
           imwrite(Image16,NameIm,'tif','Compression','none'); 
           clear NameIm Image16 

       if Opt_Tcorr
            TAv(i1)=trimmean(Ts(:),20);
       end

clear Ts 
    %  catch 
    %         display(['Issue with image ' DirIm(i1,1).name(1:end-4)]);
    %  end   
    end


    %Temperature correction, if needed
    if Opt_Tcorr

        %% SmoothingFunction 
        % From Trimmean, calculate mean temperature 
        T_Smooth = nan(size(TAv));
        for i1=1:length(T_Smooth)
            T_Smooth(i1) = nanmean(TAv(max(1,i1-Par.SmoothLength*0.5):min(length(T_Smooth),i1+Par.SmoothLength*0.5)));
        end
        %Correction
        TCorr = nanmean(T_Smooth)-T_Smooth;

        %% Create directory
        LocCorr = [LocExp 'Ts_Corr\']; mkdir(LocCorr);

        %% loop: load photos, and perform the correction
        for i1=1:length(DirIm)
            %Load images
            A = imread([LocExp DirIm(i1,1).name]);A=double(A);
            %Convert back to Ts
            Ts= (Par.TLimHigh-Par.TLimLow).*A./65535 + Par.TLimLow;
            %Convert images to 16-bit tif
            TsConv = (Ts + TCorr(i1) -Par.TLimLow)./(Par.TLimHigh-Par.TLimLow).*65535;
            Image16 = uint16(round(TsConv));
               NameIm = [LocCorr DirIm(i1,1).name];
               imwrite(Image16,NameIm,'tif','Compression','none');
        end
    end

    %Output metadata
    OutText = {'Tlow_Tif'; 'Thigh_Tif'; 'emissivity'; 'TBg'; 'Par.Tair'; 'RelHum';'Par.Tau_atm'; 'FlightHeight'; 'Opt_Sharp'; 'Opt_Tcorr'; 'Opt_Tbr'};
    OutDat = {Par.TLimLow;Par.TLimHigh;Par.emiss;Par.Tbg;Par.Tair;Par.RH;Par.Tau_atm;Par.FH;Opt_Sharp; Opt_Tcorr; Opt_Tbr};
    OutT = table([OutText OutDat]);
    writetable(OutT,[LocExp 'Parameters.xlsx']);
    if Opt_Tcorr
        OutText = {'Tlow_Tif'; 'Thigh_Tif'; 'emissivity'; 'TBg'; 'Par.Tair'; 'RelHum';'Par.Tau_atm'; 'FlightHeight'; 'Opt_Sharp'; 'Opt_Tcorr'; 'Opt_Tbr'; 'Par.SmoothLength'};
        OutDat = {Par.TLimLow;Par.TLimHigh;Par.emiss;Par.Tbg;Par.Tair;Par.RH;Par.Tau_atm;Par.FH;Opt_Sharp; Opt_Tcorr; Opt_Tbr; Par.SmoothLength};
        writetable(OutT,[LocCorr 'Parameters.xlsx']);
    end
    %% Sharpening data

    if Opt_Sharp
        LocUnsharp = [LocExp 'Unsharp\'];mkdir(LocUnsharp);
        if Opt_Tcorr
            LocUnsharpCorr = [LocCorr 'Unsharp\'];mkdir(LocUnsharpCorr);
        end

        % Step 1: Define sharpness of all images

     % DirIm = dir([LocExp '*.tif']);
      if Opt_Csv
          NameTif = struct2cell(DirIm)';NameTif = NameTif(:,1);
          for i3 =1:length(NameTif)
              NameTif{i3,1}=[NameTif{i3,1}(1:end-4) '.tif'];
          end
      else
        NameTif = struct2cell(DirIm)';NameTif = NameTif(:,1);
      end
% 
%         TSharp = nan(length(DirIm), 1);Gaps = TSharp;
% 
%         for i1=1:length(TSharp)
%             %Open image
%             A =imread([LocExp DirIm(i1,1).name(1:end-4) '.tif']);
%             A = double(A);
%             [TSharp(i1)]=estimate_sharpness(A);
%         end
%         clear i1 A

        %% Step 2 (sharpening): Define optimal sharpness

        Par.PTileThresh_Sel=Par.PTileThresh;
        x=0;
        while x==0
                Thresh = prctile(TSharp,Par.PTileThresh_Sel);
                [~,IDEdge] = min(abs(TSharp-Thresh));
                A =imread([LocIm DirIm(IDEdge,1).name(1:end-4) '.tif']);
                A=double(A);LL = min(A(:)); ML =max(A(:));
                Handles.H=figure('Name','Sharpness');  imshow(A,[LL ML]);
                prompt = 'Sharp enough? (1/0)';
                x = input(prompt);
                close(Handles.H);
            if x==0
                Par.PTileThresh_Sel = Par.PTileThresh_Sel + Par.StepThresh;
            end
        end

        ImageDropped = sum(TSharp < Thresh);
        display(['Selected threshold at percentile ' num2str(Par.PTileThresh_Sel) ' (at sharpness value of ' num2str(Thresh) '); ' num2str(ImageDropped) ' images of total of ' num2str(length(TSharp)) ' dropped']);

        %% Step 3 (sharpening): select images to move
        %Calculate GAPS
        Gaps = zeros(size(TSharp));
        for i1=ceil(Par.MaxGap*0.5):length(DirIm)-ceil(Par.MaxGap*0.5)
            if TSharp(i1) < Thresh
                R = 1;
                while (i1+R-1 < length(DirIm)-ceil(Par.MaxGap*0.5)) && TSharp(i1+R-1) < Thresh
                    R=R+1;
                end
                Gaps(i1)=R;
                clear R
            end
        end

        MoveFile = zeros(size(TSharp));
        MoveFile(TSharp < Thresh)=1;
        for i1=ceil(Par.MaxGap*0.5):length(DirIm)-ceil(Par.MaxGap*0.5)
            if Gaps(i1) > Par.MaxGap & (TSharp(i1)==max(TSharp(max(1,i1-floor(Par.MaxGap*0.5)):min(length(DirIm),i1+floor(Par.MaxGap*0.5)))))
                MoveFile(i1)=0;
            end
        end


        %% Step 4 (sharpening): move tifs
        MoveTif = NameTif(MoveFile==1);
        if ~isempty(MoveTif)
            for i2=1:length(MoveTif)
                a=movefile([LocExp MoveTif{i2}],[LocUnsharp MoveTif{i2}]);
                if Opt_Tcorr
                    a=movefile([LocCorr MoveTif{i2}],[LocUnsharpCorr MoveTif{i2}]);
                end
            end
        end
    end
end

%% Part on estimating camera positions based on multispectral data

if ~NoGPS

    TTime = nan(length(DirIm),1);

    % Look for the meta-file 

    % DirMeta = dir([LocIm '*_meta.csv']);
    % if round(length(DirIm)/10,0) == round(length(DirMeta)/10,0) %If by accident you saved a metadata file for each image, this part corrects for this.
    %    for k1=1:length(DirMeta)
    %        if k1==1 
    %            ThermalInf = ImportMeta([LocIm DirMeta(1,1).name]);
    %        else
    %            ThermalInf = [ThermalInf; ImportMeta([LocIm DirMeta(k1,1).name])];
    %        end
    %    end
    % else
    %     ThermalInf = ImportMeta([LocIm DirMeta(1,1).name]);
    % end
    if Opt_PosMulti

        TimeDiff = 0;
        Bias=[0,0, 0,0,0,0]; 
        %try
          
            if length(DirIm) == length(ThermalInf.Date)
            

                %% Get the name, time & original GPS of the multispectral data files (MultiFile)
                Add= [LocMulti '/**/*_1.tif']; %Scan in folders and subfolders
                ListMulti = dir(Add);
                 if isempty(ListMulti)
                     Add = [LocMulti '/**/*_1.tif'];
                     ListMulti = dir(Add);clear Add
                 end
                  ListMulti=struct2cell(ListMulti)';ListMulti=ListMulti(:,1:2); %date and time on Listmulti not always reliable, better extract from EXIF
                  MultiFile.TimeMulti = nan(size(ListMulti));
                  MultiFile.Lat= nan(size(ListMulti));%This is the original latitude of the camera - needed to identify the cameras, in case several cameras have the same name.
                  MultiFile.Lon= nan(size(ListMulti));%This is the original longitude of the camera - needed to identify the cameras, in case several cameras have the same name.
                  MultiFile.Alt= nan(size(ListMulti));%This is the original altitude of the camera - needed to identify the cameras, in case several cameras have the same name.
                  MultiFile.LabelMulti = cell(size(ListMulti));
                  MultiFile.ListMulti = ListMulti; clear ListMulti % MultiFile refers to the multispectral camera list (in the file)
                  for i2=1:length(MultiFile.ListMulti)
                                ImMulti = imfinfo([MultiFile.ListMulti{i2,2} '/' MultiFile.ListMulti{i2,1}]); DT = ImMulti.DateTime;
                                MultiFile.TimeMulti(i2) = datenum(str2num(DT(1:4)),str2num(DT(6:7)),str2num(DT(9:10)),...
                                                        str2num(DT(12:13)),str2num(DT(15:16)),str2num(DT(18:19)));
                               MultiFile.LabelMulti{i2} = str2num(MultiFile.ListMulti{i2}(5:8));
                        %Read latitude information from camera
                               if strcmp(ImMulti.GPSInfo.GPSLatitudeRef, 'S')
                                HLat = -1*(ImMulti.GPSInfo.GPSLatitude(1) + ImMulti.GPSInfo.GPSLatitude(2)/60 + ImMulti.GPSInfo.GPSLatitude(3)/3600);
                               else
                                   HLat = ImMulti.GPSInfo.GPSLatitude(1) + ImMulti.GPSInfo.GPSLatitude(2)/60 + ImMulti.GPSInfo.GPSLatitude(3)/3600;
                               end
                               MultiFile.Lat(i2)=HLat; clear HLat
                        %Read longitude information from camera
                               if strcmp(ImMulti.GPSInfo.GPSLongitudeRef, 'W')
                                HLon = -1*(ImMulti.GPSInfo.GPSLongitude(1) + ImMulti.GPSInfo.GPSLongitude(2)/60 + ImMulti.GPSInfo.GPSLongitude(3)/3600);
                               else
                                   HLon = ImMulti.GPSInfo.GPSLongitude(1) + ImMulti.GPSInfo.GPSLongitude(2)/60 + ImMulti.GPSInfo.GPSLongitude(3)/3600;
                               end
                               MultiFile.Lon(i2)=HLon; clear HLon
                         %Read Altitude
                               MultiFile.Alt(i2)=ImMulti.GPSInfo.GPSAltitude(1); 

                  end
                  clear Add 

                  %% Load the camera positions as exported from Metashape (MultiPos)

                      MultiPos = ImportMulti_Pos([LocTxtMulti NameTxt]);

%                  % MultiPos.Names = NameMulti_Pos.MultiPos; %clear NameMulti_Pos
%                   MultiPos.X_est= MultiPos.X_est.*10^-6;
%                   MultiPos.Y_est= MultiPos.Y_est.*10^-6;
%                   MultiPos.Z_est= MultiPos.Z_est.*10^-6;
%                   MultiPos.Yaw_est= MultiPos.Yaw_est.*10^-6;
%                   MultiPos.Pitch_est= MultiPos.Pitch_est.*10^-6;
%                   MultiPos.Roll_est= MultiPos.Roll_est.*10^-6;


               %   NameMulti_Pos = char(NameMulti_Pos);

           %% Match both multi lists: in practice, create a list with times, and corresponding (corrected) locations/pitch/yaw/roll data
                  % A simple ismember function no longer works, since there
                  % are several cameras with the same name. We therefore
                  % identify each camera separately, based on the original
                  % location.
                  MultiPos.Time = nan(size(MultiPos.XLongitude));  
                  for i1=1:length(MultiPos.XLongitude)
                      IdMF = find(strcmp(MultiFile.ListMulti(:,1),string(MultiPos.Label(i1))));
                      %Information of the GPS location apparently doesn't
                      %work. We'll have to use the thermal camera data: 
                      %1) For the cameras identified, find the thermal
                      %camera with the same time, and record its location
                      if isempty(IdMF)
                          display(['Issue with finding data on Camera ' string(MultiPos.Label(i1))]);
                      elseif length(IdMF) ==1 %only one camera, so you can immediately notate the right data for this camera
                            MultiPos.Time(i1) = MultiFile.TimeMulti(IdMF)
                      else %Several cameras with the same name...we need to find out which one couples to which
%                         Dist = nan(length(IdMF),1);
%                         Distb= nan(length(IdMF),1);
                        Distc= nan(length(IdMF),1);
                        for i2=1:length(IdMF)
%                                Dist(i2) = acos(cosd(90-MultiPos.YLatitude(i1)) .* cosd(90-MultiFile.Lat(IdMF(i2))) + sind(90-MultiPos.YLatitude(i1)) .* sind(90-MultiFile.Lat(IdMF(i2)))) .* cosd(MultiPos.XLongitude(i1)-MultiFile.Lon(IdMF(i2))) * 3958.76;
%                                %From https://nl.mathworks.com/matlabcentral/answers/435843-calculate-distance-from-gps-coordinates
%                                [Distb(i2) ~]=lldistkm([MultiPos.YLatitude(i1) MultiPos.XLongitude(i1)],[MultiFile.Lat(IdMF(i2)) MultiFile.Lon(IdMF(i2))]);
                               Distc(i2) = stdist([MultiPos.YLatitude(i1) MultiFile.Lat(IdMF(i2))],[MultiPos.XLongitude(i1) MultiFile.Lon(IdMF(i2))]);%In degrees, see https://nl.mathworks.com/help/map/ref/stdist.html 
                        end
                        [~,MD] = min(Distc);
                        MultiPos.Time(i1) = MultiFile.TimeMulti(IdMF(MD));
                      end                     
                  end
                  clear MD Distc i2 IdMF
                %Finally, sort according to MultiPos.Time
                [~,a] = sort(MultiPos.Time);
                MultiSorted.Time = MultiPos.Time(a);
                MultiSorted.Label = MultiPos.Label(a);
                MultiSorted.X_est = MultiPos.X_est(a);MultiSorted.Y_est = MultiPos.Y_est(a);MultiSorted.Z_est = MultiPos.Z_est(a);
                MultiSorted.Yaw_est = MultiPos.Yaw_est(a);MultiSorted.Pitch_est = MultiPos.Pitch_est(a);MultiSorted.Roll_est = MultiPos.Roll_est(a);
                
                %Remove Non-Nan values
                ID = ~isnan(MultiPos.X_est);
               % MultiSorted = MultiSorted(ID==1,:);
                MultiSorted.Time=MultiSorted.Time(ID);MultiSorted.Label =MultiSorted.Label(ID); MultiSorted.X_est =MultiSorted.X_est(ID);
                MultiSorted.Y_est =MultiSorted.Y_est(ID);MultiSorted.Z_est =MultiSorted.Z_est(ID);MultiSorted.Yaw_est =MultiSorted.Yaw_est(ID);
                MultiSorted.Roll_est =MultiSorted.Roll_est(ID);MultiSorted.Pitch_est =MultiSorted.Pitch_est(ID);
                clear ID
                
                %No 2 images at the same second
                MultiSorted.Time_corr =  MultiSorted.Time;
                    for i6=2:length(MultiSorted.Time)
                        if MultiSorted.Time_corr(i6)==MultiSorted.Time_corr(i6-1)
                            MultiSorted.Time_corr(i6)=MultiSorted.Time_corr(i6)+0.5./(24.*3600);
                        end
                    end

                  %% Match the thermal data with the multispectral position, based on the time stamp on the thermal data
                        NameThermal = struct2cell(DirIm);NameThermal=NameThermal';NameThermal = NameThermal(:,1);
                        X = nan(length(NameThermal),1);Y=X;Z=X;yaw = X; pitch=X; roll=X;
                     TableThermal=table(NameThermal,X,Y,Z,yaw,pitch,roll);clear X Y pitch yaw Z roll
                     for i4=1:length(TableThermal.X)
                         try
                            TableThermal.X(i4) = interp1(MultiSorted.Time_corr, MultiSorted.X_est, TTime_corr(i4)-TimeDiff/(24*3600)) + Bias(1);
                            TableThermal.Y(i4) = interp1(MultiSorted.Time_corr(:), MultiSorted.Y_est, TTime_corr(i4)-TimeDiff/(24*3600))+ Bias(2);
                            TableThermal.Z(i4) = interp1(MultiSorted.Time_corr(:), MultiSorted.Z_est, TTime_corr(i4)-TimeDiff/(24*3600))+ Bias(3);
                            TableThermal.yaw(i4) = interp1(MultiSorted.Time_corr, MultiSorted.Yaw_est, TTime_corr(i4)-TimeDiff/(24*3600))+ Bias(4);
                            TableThermal.pitch(i4) = interp1(MultiSorted.Time_corr, MultiSorted.Pitch_est, TTime_corr(i4)-TimeDiff/(24*3600))+ Bias(5);
                            TableThermal.roll(i4) = interp1(MultiSorted.Time_corr, MultiSorted.Roll_est, TTime_corr(i4)-TimeDiff/(24*3600))+ Bias(6);
                         end
                     end

                     %Look for NaN-values, and replace them with mean value
                     %from above or below
                     A = find(isnan(TableThermal.X));
                     if ~isempty(A)
                         for k1=1:length(A)  
                            try
                            TableThermal.X(A(k1)) = nanmean(TableThermal.X(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));
                            TableThermal.Y(A(k1)) = nanmean(TableThermal.Y(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));
                            TableThermal.Z(A(k1)) = nanmean(TableThermal.Z(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));
                            TableThermal.yaw(A(k1)) = nanmean(TableThermal.yaw(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));
                            TableThermal.pitch(A(k1)) = nanmean(TableThermal.pitch(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));
                            TableThermal.roll(A(k1)) = nanmean(TableThermal.roll(max(1,A(k1)-5):min(A(k1)+5, length(TableThermal.X))));

                        end
                     end


              Handles.H=figure('Name','XY-pattern Thermal');  scatter(TableThermal.X,TableThermal.Y, 'o');
                prompt = 'Correct XY-values? (based on multispectral data) (1/0)';
                x = input(prompt);
                close(Handles.H);
            if x==1
                writetable(TableThermal,[LocExp 'ThermalPos_FromMulti.txt']);
            else
                display('Issue with matching thermal and multispectral image position: Data not saved');
            end


            else
                display('Issue with matching thermal and multispectral image position: metadata thermal not correct');
                     end
            end
%          catch
%               display('Issue with matching thermal and multispectral image position')
%          end
    end

%if no data from multispectral flights are available, get XYZ from the GNSS of the thermal camera.
    if length(DirIm) == length(ThermalInf.Date)
            NameThermal = struct2cell(DirIm);NameThermal=NameThermal';NameThermal = NameThermal(:,1);
            
            TableThermal=table(NameThermal,ThermalInf.Latitude,ThermalInf.Longitude,ThermalInf.Altitude,ThermalInf.Yaw,ThermalInf.Pitch,ThermalInf.Roll, 'VariableNames', {'Names', 'lat','lon','alt', 'yaw','pitch','roll',});
                        
            Handles.H=figure('Name','XY-pattern Thermal');  scatter(TableThermal.lat,TableThermal.lon, 'o');
            prompt = 'Correct XY-values? (based on thermal GNSS sensor) (1/0)';
            x = input(prompt);
            close(Handles.H);
            writetable(TableThermal,[LocExp 'ThermalPos_FromThermal.txt']);
    else
       display('Issue with GNSS position thermal camera: metadata thermal not correct');
    end

end


%% 

function [sharpness]=estimate_sharpness(G)

[Gx, Gy]=gradient(G);

S=sqrt(Gx.*Gx+Gy.*Gy);
sharpness=sum(sum(S))./(numel(Gx));

end

%%
function A= f_readIm(ImName, Opt_CSV)

  %Open data
   if Opt_CSV
     A=readmatrix(ImName);
   else
    A = imread(ImName);A=double(A); A=A.*0.04-273.15;
   end
end

% 
% function namemeta = ImportMeta(filename)
% %IMPORTFILE Import data from a text file
% %  NAMEMETA = IMPORTFILE(FILENAME) reads data from text file FILENAME
% %  for the default selection.  Returns the data as a table.
% %
% %  NAMEMETA = IMPORTFILE(FILE, DATALINES) reads data for the specified
% %  row interval(s) of text file FILENAME. Specify DATALINES as a
% %  positive scalar integer or a N-by-2 array of positive scalar integers
% %  for dis-contiguous row intervals.
% %
% %  Example:
% %  namemeta = importfile("E:\2021_drone\20210814_Gedinne\Gedinne1\thermal\testje\name_meta.csv", [2, Inf]);
% %
% %  See also READTABLE.
% %
% % Auto-generated by MATLAB on 25-Nov-2021 12:04:15
% 
% %% Input handling
% 
% % If dataLines is not specified, define defaults
% %if nargin < 2
%     dataLines = [2, Inf];
% %end
% 
% %% Setup the Import Options and import the data
% opts = delimitedTextImportOptions("NumVariables", 14);
% 
% % Specify range and delimiter
% opts.DataLines = dataLines;
% opts.Delimiter = ";";
% 
% % Specify column names and types
% opts.VariableNames = ["Date", "Time", "Latitude", "Longitude", "Pitch", "Roll", "Yaw", "Altitude", "Satellites", "TrackHeading", "MagneticHeading", "TempHousing", "TempFPA", "FFC_AfterFrame"];
% opts.VariableTypes = ["categorical", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Specify variable properties
% opts = setvaropts(opts, "Date", "EmptyFieldRule", "auto");
% opts = setvaropts(opts, "Time", "InputFormat", "HH:mm:ss");
% 
% % Import the data
% namemeta = readtable(filename, opts);
% 
% end


function namemeta = ImportMeta(filename)
%IMPORTFILE Import data from a text file
%  NAMEMETA = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  NAMEMETA = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  namemeta = importfile("E:\2021_drone\20210814_Gedinne\Gedinne1\thermal\testje\name_meta.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 25-Nov-2021 12:04:15

%% Input handling

% If dataLines is not specified, define defaults
%if nargin < 2
    dataLines = [2, Inf];
%end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 14);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["Date", "Time", "Latitude", "Longitude", "Pitch", "Roll", "Yaw", "Altitude", "Satellites", "TrackHeading", "MagneticHeading", "TempHousing", "TempFPA", "FFC_AfterFrame"];
opts.VariableTypes = ["categorical", "datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Date", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "Time", "InputFormat", "HH:mm:ss");

% Import the data
namemeta = readtable(filename, opts);

end

function NameMulti_Pos =Imp_MultiNames(filename)
%IMPORTFILE Import data from a text file
%  GED1MULTIPOS = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a cell array.
%
%  GED1MULTIPOS = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  Ged1MultiPos = importfile("E:\2021_drone\20210814_Gedinne\Agisoft\210814_Ged1_Multi_Pos.txt", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 25-Nov-2021 12:51:37

%% Input handling

% If dataLines is not specified, define defaults

    dataLines = [2, Inf];


%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Label", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13"];
opts.SelectedVariableNames = "Label";
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Label", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Label", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13"], "EmptyFieldRule", "auto");

% Import the data
Ged1MultiPos = readtable(filename, opts);

%% Convert to output type
NameMulti_Pos = table2cell(Ged1MultiPos);
% numIdx = cellfun(@(x) ~isnan(str2double(x)), Ged1MultiPos);
% NameMulti_Pos(numIdx) = cellfun(@(x) {str2double(x)}, Ged1MultiPos(numIdx));
end




function MultiPos = ImportMulti_Pos(filename, dataLines)


%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [3, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Label", "XLongitude", "YLatitude", "ZAltitude", "Yaw", "Pitch", "Roll", "X_est", "Y_est", "Z_est", "Yaw_est", "Pitch_est", "Roll_est"];
opts.VariableTypes = ["categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Label", "EmptyFieldRule", "auto");

% Import the data
MultiPos = readtable(filename, opts);
IDKeep = ~isnan(MultiPos.XLongitude);
MultiPos = MultiPos(IDKeep,:);

end

function [d1km d2km]=lldistkm(latlon1,latlon2)
% format: [d1km d2km]=lldistkm(latlon1,latlon2)
% Distance:
% d1km: distance in km based on Haversine formula
% (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)
% d2km: distance in km based on Pythagoras’ theorem
% (see: http://en.wikipedia.org/wiki/Pythagorean_theorem)
% After:
% http://www.movable-type.co.uk/scripts/latlong.html
%
% --Inputs:
%   latlon1: latlon of origin point [lat lon]
%   latlon2: latlon of destination point [lat lon]
%
% --Outputs:
%   d1km: distance calculated by Haversine formula
%   d2km: distance calculated based on Pythagoran theorem
%
% --Example 1, short distance:
%   latlon1=[-43 172];
%   latlon2=[-44  171];
%   [d1km d2km]=distance(latlon1,latlon2)
%   d1km =
%           137.365669065197 (km)
%   d2km =
%           137.368179013869 (km)
%   %d1km approximately equal to d2km
%
% --Example 2, longer distance:
%   latlon1=[-43 172];
%   latlon2=[20  -108];
%   [d1km d2km]=distance(latlon1,latlon2)
%   d1km =
%           10734.8931427602 (km)
%   d2km =
%           31303.4535270825 (km)
%   d1km is significantly different from d2km (d2km is not able to work
%   for longer distances).
%
% First version: 15 Jan 2012
% Updated: 17 June 2012
%--------------------------------------------------------------------------

radius=6371;
lat1=latlon1(1)*pi/180;
lat2=latlon2(1)*pi/180;
lon1=latlon1(2)*pi/180;
lon2=latlon2(2)*pi/180;
deltaLat=lat2-lat1;
deltaLon=lon2-lon1;
a=sin((deltaLat)/2)^2 + cos(lat1)*cos(lat2) * sin(deltaLon/2)^2;
c=2*atan2(sqrt(a),sqrt(1-a));
d1km=radius*c;    %Haversine distance

x=deltaLon*cos((lat1+lat2)/2);
y=deltaLat;
d2km=radius*sqrt(x*x + y*y); %Pythagoran distance

end