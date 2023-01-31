%{
Online supplementary materials of the paper titled 
"Distributionally Robust State Estimation for Jump Linear Systems"
Authored By: Shixiong Wang
From the Institute of Data Science, National University of Singapore

@Author: Shixiong Wang
@Date: 16 May 2022
@Email: s.wang@u.nus.edu; wsx.gugo@gmail.com
@Site: https://github.com/Spratm-Asleaf/DRSE-Jump
%}

%{
    Parse GNSS data provided by the UniStrong Co., Ltd., Beijing 
    The "FVI" statement is used
%}

clear all;
close all;
clc;

%% Open source files
fidGPS = fopen('GPSRaw.dat', 'r');             % Concentional GPS data
if fidGPS <= 0
    disp('Error in opening file: GPS');
    return;
end

fidRTK = fopen('RTKRaw.dat', 'r');             % Very-high-accuracy GPS data
if fidRTK <= 0
    disp('Error in opening file: RTK');
    return;
end

%% Retrived Data: For the same trajectory
% Records by Conventional GPS
LLH_GPS = [];           % Longitude-Latitude-Height Coordinate
XYZ_GPS = [];           % X, Y, Z in World Geodetic System (WGS) Coordinate; (see: https://en.wikipedia.org/wiki/World_Geodetic_System)
ENU_GPS = [];           % X, Y, Z in Local Tangent Plane Coordinate; (see: https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates)
TIME_GPS = [];          % Time
GPSStatus = [];         % GPS Status: 1 - Normal GPS Working Mode
% Records by RTK
LLH_RTK = [];
XYZ_RTK = [];
ENU_RTK = [];
TIME_RTK = [];          % Time
VEL_RTK = [];           % Real Moving Velocities Returned by RTK
RTKStatus = [];         % GPS Status: 4 - Reliable RTK Working Mode

%% Origin of the ENU coordinate
% The first point from RTK will be treated as the origin
phi0 = 0;               % Longitude
lamda0 = 0;             % Latitude

%% Parse RTK
k = 0;
while ~feof(fidRTK)
    msg = fgetl(fidRTK);
    msg = [msg ','];
    len = length(msg);
    
    if len <= 10
        continue;
    else
        if strcmp(msg(1:10), '$PSAT,FVI,') == 0
            continue;
        else
            k = k + 1;
        end
    end

    data = cell(31,1);

    substr = [];

    j = 1;
    for i = 1:len
        if msg(i) == '*'
            msg(i) = ',';
        end
        
        if msg(i) ~= ','
            substr = [substr msg(i)];
            continue;
        end

        if strcmp(substr, '$PSAT') == 1
            data{j} = substr;
        elseif strcmp(substr, 'FVI') == 1
            data{j} = [data{j} '-FVI'];
            j = j+1;
        else
            data{j} = substr;
            j = j+1;
        end
        
        substr = [];
    end
    
    %% Status
    RTKStatus = [RTKStatus; str2double(data(27))];
    
    %% Time
    TIME_RTK = [TIME_RTK str2double(data{2})];

    %% LLH: Long-Lat-Height
    phi = str2double(data{3})*pi/180;   % longitude
    lamda = str2double(data{4})*pi/180; % latitude
    height = str2double(data{5});
    
    LLH_RTK = [LLH_RTK [phi; lamda; height]];

    %% XYZ: GWS84
    a = 6378137.00;
    e = 0.081819191; %f = 1/298.257223563;
    R = a./sqrt(1-e^2 * sin(phi)^2);

    XYZ_RTK = [XYZ_RTK [
                (R + height)*cos(phi)*cos(lamda);
                (R + height)*cos(phi)*sin(lamda);
                (R * (1-e^2) + height)*sin(phi);
                ];
    ];

    %% ENU
    if k == 1
        phi0 = phi;
        lamda0 = lamda;
        
        ENU_RTK = [ENU_RTK [0;0;0]];
        VEL_RTK = [VEL_RTK [0;0;0]];
        
        continue;
    else
        M = [
            -sin(lamda0)                cos(lamda0)             0
            -sin(phi0)*cos(lamda0)     -sin(phi0)*sin(lamda0)   cos(phi0)
             cos(phi0)*cos(lamda0)      cos(phi0)*sin(lamda0)   sin(phi0)
        ];
        enu = M*(XYZ_RTK(:,k) - XYZ_RTK(:,1));
        
        ENU_RTK = [ENU_RTK enu];
    end
    
    %% Real Velocity
    vel = [str2double(data{15});str2double(data{16});str2double(data{17})];
    VEL_RTK = [VEL_RTK vel];
end

%% Parse GPS
k = 0;
while ~feof(fidGPS)
    msg = fgetl(fidGPS);
    msg = [msg ','];
    len = length(msg);
    
    if len <= 10
        continue;
    else
        if strcmp(msg(1:10), '$PSAT,FVI,') == 0
            continue;
        else
            k = k + 1;
        end
    end

    data = cell(31,1);

    substr = [];

    j = 1;
    for i = 1:len
        if msg(i) == '*'
            msg(i) = ',';
        end
        
        if msg(i) ~= ','
            substr = [substr msg(i)];
            continue;
        end

        if strcmp(substr, '$PSAT') == 1
            data{j} = substr;
        elseif strcmp(substr, 'FVI') == 1
            data{j} = [data{j} '-FVI'];
            j = j+1;
        else
            data{j} = substr;
            j = j+1;
        end
        
        substr = [];
    end
    
    %% Status
    GPSStatus = [GPSStatus; str2double(data(27))];
    
    %% Time
    TIME_GPS = [TIME_GPS str2double(data{2})];

    %% LLH: Long-Lat-Height
    phi = str2double(data{3})*pi/180;   % long
    lamda = str2double(data{4})*pi/180; % lat
    height = str2double(data{5});
    
    LLH_GPS = [LLH_GPS [phi; lamda; height]];

    %% XYZ: GWS84
    a = 6378137.00;
    e = 0.081819191; %f = 1/298.257223563;
    R = a./sqrt(1-e^2 * sin(phi)^2);

    XYZ_GPS = [XYZ_GPS [
                (R + height)*cos(phi)*cos(lamda);
                (R + height)*cos(phi)*sin(lamda);
                (R * (1-e^2) + height)*sin(phi);
                ];
    ];

    %% ENU
    M = [
        -sin(lamda0)                cos(lamda0)             0
        -sin(phi0)*cos(lamda0)     -sin(phi0)*sin(lamda0)   cos(phi0)
         cos(phi0)*cos(lamda0)      cos(phi0)*sin(lamda0)   sin(phi0)
    ];
    enu = M*(XYZ_GPS(:,k) - XYZ_RTK(:,1));

    ENU_GPS = [ENU_GPS enu];
end

fclose(fidGPS);
fclose(fidRTK);

%% Write Retrived Data;
fid = fopen('..\main\Data.txt', 'w');
len = length(TIME_GPS);
fprintf(fid, 'Time\t\tGPS_X\t\tGPS_Y\t\tRTK_X\t\tRTK_Y\t\tRTK_Vx\t\tRTK_Vy\n');
for i = 1:len
    fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\n', 0.1*(i-1), ENU_GPS(1,i), ENU_GPS(2,i), ENU_RTK(1,i), ENU_RTK(2,i), VEL_RTK(1,i), VEL_RTK(2,i));
end
fclose(fid);

%% Plot RTK Trajectory (Real Trajectory)
% plot3(ENU_GPS(1,:),ENU_GPS(2,:),ENU_GPS(3,:));
% plot3(ENU_RTK(1,:),ENU_RTK(2,:),ENU_RTK(3,:));