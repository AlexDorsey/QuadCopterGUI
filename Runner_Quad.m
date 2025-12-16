% RUNNER_QUAD
% Created By Alex Dorsey, adorsey4@umbc.edu, Fall 2025
% Follows standard reset-simutation and step forward format seen in
% training enviroments
%
% This script launches an interactive quadcopter simulation with:
% - 3D visualization, camera tracking, and ground shadow
% - GUI panels for dynamics, estimation, simulation, controller, and waypoints
% - Keyboard based user commands (W A S D, Q E, Shift, Control)
% - Optional waypoint navigation in ENU frame with cruise and slow radius logic
% - Inner loop attitude control and optional adaptive DMAC controller
% - Roll and pitch estimation via Kalman filter with optional lowpass filter
% - Logging of attitude, position, and theta FFT for post processing
%
% x(1): phi: pitch
% x(2): theta: roll
% x(3): psi: yaw
% x(4:6): body rates [p;q;r]
% x(7): east
% x(8): north
% x(9): up
% x(10): v_east
% x(11): v_north
% x(12): v_up

clc; clear; close all;
addpath('..\'); rng(123);

SIMDATA = resetSim();
SIMDATA = runSimStep(SIMDATA);  
setappdata(SIMDATA.figureHandle,'SIMDATA',initUiState(SIMDATA)); 

while isvalid(SIMDATA.figureHandle)
    % read UI flags set by buttons
    ui = getappdata(SIMDATA.figureHandle,'SIMDATA');
    if isempty(ui)
        ui = initUiState(SIMDATA); 
        setappdata(SIMDATA.figureHandle,'SIMDATA',ui); 
    end

    % restart requested
    if ui.requestReset
        SIMDATA = resetSim(false);
        ui = initUiState(SIMDATA);
        ui.running = true;
        setappdata(SIMDATA.figureHandle,'SIMDATA',ui);
        continue
    end

    % run one step only when running
    if ui.running
        SIMDATA = runSimStep(SIMDATA);   
    else
        pause(0.01);                     
    end

    drawnow limitrate
end

%% Main Local Functions
function ui = initUiState(SIMDATA)
    ui.running = false;          % Start sets true
    ui.requestReset = false;     
    ui.hFig = SIMDATA.figureHandle;
end

function SIMDATA = resetSim(remakeFigureFlag)
    persistent OBJs axesHandle figureHandle
    %% User Parameters
    USER_PARAMS.userCommandedControl = true;
    USER_PARAMS.AngleCommandFilterTau = 0.25;
    USER_PARAMS.maxVelocity = 100;
    USER_PARAMS.acceleration = 100; 
    USER_PARAMS.angular_velocity = pi/2; 
    USER_PARAMS.xLim = 1000;
    % plotting
    USER_PARAMS.PlotRealTime = true;
    USER_PARAMS.logData = false;
    USER_PARAMS.DrawShadow = true;
    USER_PARAMS.fps = 30;
    USER_PARAMS.quadScale = 5;
    USER_PARAMS.propSpinRPM = 1000;
    % control system
    USER_PARAMS.UseAdaptiveController = true;
    USER_PARAMS.EnableWaypoint = false;
    USER_PARAMS.Waypoints = [10 20 20 10; 10 10 20 20; 10 20 30 40];  % 3 x N [east; north; up]
    USER_PARAMS.waypointCruiseSpeed = 5;   % m/s
    USER_PARAMS.waypointSlowRadius  = 20;    % m  
    % drone physical parameters
    USER_PARAMS.mass_kg = 1;
    USER_PARAMS.Jx = 0.5; USER_PARAMS.Jy = 0.5; USER_PARAMS.Jz = 0.25;
    
    %% Simulation
    dynamics = @(x, moments, force, params) quad_dynamics(x, moments, force, params);
    
    dt = 1/150;
    
    moment = zeros(3,1); 
    stateVector = zeros(12,1);
    stateVector(9,1) = 5;
    
    if exist('remakeFigureFlag','var') && remakeFigureFlag == false
        % do nothing
    else
        [OBJs, axesHandle, figureHandle] = initQaudFrameAndGUI(USER_PARAMS);
    end
    DMAC_STATES = initDMAC_Settings(dt); % DMAC SETTINGS, go into this function to change DMAC Settings

    if USER_PARAMS.userCommandedControl
        setupKeyHandler(figureHandle);
        userCommands.velocityX = 0;
        userCommands.velocityY = 0;
        userCommands.velocityZ = 0;
        userCommands.psi = 0;
    end
    SIMDATA.prevResponseVel = [0;0;0];
    SIMDATA.USER_PARAMS = USER_PARAMS;
    SIMDATA.OBJs = OBJs;
    SIMDATA.OBJs.AttitudeLog = [];
    SIMDATA.OBJs.PositionLog = [];
    SIMDATA.axesHandle = axesHandle;
    SIMDATA.userCommands = userCommands;
    SIMDATA.stateVector = stateVector;
    SIMDATA.moment = moment;
    SIMDATA.DMAC_STATES = DMAC_STATES;
    SIMDATA.dynamics = dynamics;
    SIMDATA.dt = dt;
    SIMDATA.figureHandle = figureHandle;
    SIMDATA.kIter = 1;
    SIMDATA.time = 0;
    processNoise = OBJs.GUI_Params.KalmanProcessNoiseField.Value;
    measurementNoise = OBJs.GUI_Params.KalmanMeasurementNoiseField.Value;
    initAngleUncertainty = OBJs.GUI_Params.KalmanInitalAngleErrorField.Value;
    KalmanFilterRollPitch(dt, [0;0;0], [0;0;9.81], true, initAngleUncertainty, measurementNoise, processNoise);
    getWaypointVelocityCommandFromWaypointList(true, zeros(3,1), zeros(12,1), USER_PARAMS);
end

function SIMDATA = runSimStep(SIMDATA)
    % Unpack
    USER_PARAMS = SIMDATA.USER_PARAMS;
    OBJs = SIMDATA.OBJs;
    axesHandle = SIMDATA.axesHandle;
    userCommands = SIMDATA.userCommands;
    stateVector = SIMDATA.stateVector;
    moment = SIMDATA.moment;
    DMAC_STATES = SIMDATA.DMAC_STATES;
    dynamics = SIMDATA.dynamics;
    kIter = SIMDATA.kIter;
    time = SIMDATA.time;
    if time == 0
        firstRunBool = true;
    else
        firstRunBool = false;
    end
    % update state     
    
    % update mass and J from GUI
    USER_PARAMS.mass_kg = OBJs.GUI_Params.massField.Value;
    USER_PARAMS.J = OBJs.GUI_Params.JTableField.Data;
    % update kalman params from GUI
    processNoise = OBJs.GUI_Params.KalmanProcessNoiseField.Value;
    measurementNoise = OBJs.GUI_Params.KalmanMeasurementNoiseField.Value;
    initAngleUncertainty = OBJs.GUI_Params.KalmanInitalAngleErrorField.Value;
    % update OuterLoop Controller from GUI
    SIMDATA.USER_PARAMS.acceleration = OBJs.GUI_Params.OuterLoopAccelerationField.Value;
    SIMDATA.USER_PARAMS.angular_velocity = OBJs.GUI_Params.OuterLoopAngularVelocityField.Value*pi/180;
    EnableAngleControlFlag = OBJs.GUI_Params.EnableAngleControlField.Value;
    AdaptiveControlFlag = OBJs.GUI_Params.EnableAdaptiveControlField.Value;
    SIMDATA.USER_PARAMS.AngleCommandFilterTau = OBJs.GUI_Params.AngleCommandFilterTauField.Value;
    
    % simulation update
    SIMDATA.dt =  1/OBJs.GUI_Params.simHzField.Value;
    
    dt = SIMDATA.dt;
    USER_PARAMS.fps = OBJs.GUI_Params.displayHzField.Value;
    
    % update time
    time = time + dt;

    % waypoint handling
    EnableWaypointFlag = OBJs.GUI_Params.EnableWaypointField.Value;
    SIMDATA.USER_PARAMS.EnableWaypoint = strcmp(EnableWaypointFlag,'On');

    % parse text into numeric 3 x N matrix
    wpString = OBJs.GUI_Params.WaypointDataField.Value;
    waypointList = str2num(wpString); 
    
    USER_PARAMS.waypointCruiseSpeed = OBJs.GUI_Params.WaypointCruiseSpeedField.Value;
    USER_PARAMS.waypointSlowRadius  = OBJs.GUI_Params.WaypointSlowRadiusField.Value;
    
    if SIMDATA.USER_PARAMS.EnableWaypoint
        % waypoint derived velocity commands
        velocityCommands = getWaypointVelocityCommandFromWaypointList(false, waypointList, stateVector, USER_PARAMS);
        [velocityCommands] = firstOrderFilter(SIMDATA.prevResponseVel, velocityCommands, dt, 0.15);
        [force, angle_command] = velocityController(kIter, stateVector, velocityCommands, dt, 0, USER_PARAMS);
        SIMDATA.prevResponseVel = velocityCommands;
    else % not waypoint
        if strcmp(EnableAngleControlFlag,'Off')
            [userCommands, velocityCommands, psi_cmd] = getKeysAndVelocityCommands(userCommands, stateVector, USER_PARAMS, dt);
            [force, angle_command] = velocityController(kIter, stateVector, velocityCommands, dt, psi_cmd, USER_PARAMS); 
        else
            [userCommands, angle_command, force] = getKeysAndAngleCommands(userCommands, stateVector, USER_PARAMS, dt);
        end
    end
    PHI = [stateVector(1:6); moment]; %[x;u]
    
    % plant propagation
    stateVector = RK4(dynamics, dt, stateVector, moment, force, USER_PARAMS);
    % sensor readings
    [~, acceleration_body, angularVelocity_body]= dynamics(stateVector, moment, force, USER_PARAMS);

    AccelWhiteNoiseMag = OBJs.GUI_Params.AccelWhiteNoiseMagField.Value;
    GyroWhiteNoiseMag = OBJs.GUI_Params.GyroWhiteNoiseMagField.Value;
    AccelBias = OBJs.GUI_Params.AccelBiasField.Value;
    GryoBias = OBJs.GUI_Params.GryoBiasField.Value;

    accel_xyz = acceleration_body + deg2rad(AccelWhiteNoiseMag)*randn(3,1) + deg2rad(AccelBias); %noise parameters
    gyro_pqr =  angularVelocity_body + deg2rad(GyroWhiteNoiseMag)*randn(3,1) + deg2rad(GryoBias); %noise parameters
    % Kalman Filtering
    
    [phi_hat, theta_hat] = KalmanFilterRollPitch(dt, gyro_pqr, accel_xyz, false, initAngleUncertainty, measurementNoise, processNoise);

    %% Inject Noise
    NoiseMagByFreq= str2num(OBJs.GUI_Params.NoiseMagFreqAngleField.Value);
    Mag = NoiseMagByFreq(1,:);
    Freq = NoiseMagByFreq(2,:);
    
    noiseSum = 0;
    for noiseIndex = 1:length(Mag)
        noiseSum = noiseSum + Mag(noiseIndex)*cos(Freq(noiseIndex)*time);
    end
    noiseSum  = noiseSum/length(Mag);
    
    phi_hat = phi_hat + noiseSum;
    theta_hat = theta_hat + noiseSum;
    
    %% lowpass filter noise
    LowpassFlag = OBJs.GUI_Params.EnableAngleLowpassField.Value;
    cutOffFreqHz = OBJs.GUI_Params.AngleLowpassCutoffField.Value;

    if strcmp(LowpassFlag,'On')
        [phi_hat, theta_hat] = lowpassFilter(firstRunBool, phi_hat, theta_hat, cutOffFreqHz, dt);
    end
    y_DMAC = stateVector(1:6);
    y_DMAC(1) = phi_hat;
    y_DMAC(2) = theta_hat;
    %% FFT
    phi_true   = stateVector(1);
    theta_true = stateVector(2);

    OBJs.AttitudeLog = [OBJs.AttitudeLog, [phi_true; phi_hat; theta_true; theta_hat]]; % LOGGING FOR PLOTTING AND fft
    
    east  = stateVector(7);
    north = stateVector(8);
    OBJs.PositionLog = [OBJs.PositionLog, [east; north]];
   

    NFFT = 512;
    if kIter > NFFT
        % phiWindow = OBJs.AttitudeLog(2,:);
        thetaWindow = OBJs.AttitudeLog(4,:);

        % phiWindow = phiWindow - mean(phiWindow);
        thetaWindow = thetaWindow - mean(thetaWindow);

        % fftPhi = fft(phiWindow);
        fftTheta = fft(thetaWindow);

        % Single sided spectrum from matlab website
        % P2Phi   = abs(fftPhi   / NFFT);
        P2Theta = abs(fftTheta / NFFT);
        % P1Phi   = P2Phi(1:NFFT/2+1);
        P1Theta = P2Theta(1:NFFT/2+1);
        Fs = 1/dt;
        fAxis = (0:(NFFT/2)) * (Fs / NFFT);

        OBJs.fAxis = fAxis;
        OBJs.P1Theta = P1Theta;
    end


    
    %% adaptive controller
    USER_PARAMS.UseAdaptiveController = strcmp(AdaptiveControlFlag,'On');
    if USER_PARAMS.UseAdaptiveController
        DMAC_STATES.Q = diag(str2num(OBJs.GUI_Params.DMAC_Qdiag_Field.Value));
        DMAC_STATES.R = diag(str2num(OBJs.GUI_Params.DMAC_Rdiag_Field.Value));
        DMAC_STATES.R0 = OBJs.GUI_Params.DMAC_R0_Field.Value;
    
        halfLifeGUI = OBJs.GUI_Params.DMAC_lambda_Field.Value;
        DMAC_STATES.lambda = 0.5^(dt / halfLifeGUI);

        [DMAC_STATES, moment] = DMAC_V2(DMAC_STATES, PHI, y_DMAC, angle_command);
        moment = clip(moment, -150,150);
        moment = moment + 1e-2*randn(3,1);
    else
        moment = innerLoopController(kIter, dt, stateVector, angle_command, USER_PARAMS);
        moment = clip(moment, -150,150);
    end
     
    % Update quad copter plot
    if USER_PARAMS.PlotRealTime
        if mod(kIter,ceil(1/(USER_PARAMS.fps*dt))) == 0 || kIter == 1
            OBJs = updateQuadPlot(kIter, axesHandle, stateVector, dt, OBJs);
            drawnow
        end
    end
        
    % if stateVector(9) < 0
    %     error("Quad hit the ground!")
    % end
        
    % repack
    SIMDATA.USER_PARAMS = USER_PARAMS;
    SIMDATA.OBJs = OBJs;
    SIMDATA.axesHandle = axesHandle;
    SIMDATA.userCommands = userCommands;
    SIMDATA.stateVector = stateVector;
    SIMDATA.moment = moment;
    SIMDATA.DMAC_STATES = DMAC_STATES;
    SIMDATA.dt = dt;
    SIMDATA.kIter = kIter + 1;
    SIMDATA.time = time;
end

function [OBJs, axesHandle, figureHandle]  = initQaudFrameAndGUI(USER_PARAMS)
    %monstrosity of a function that grew as I added functionality
    %iteratively to the GUI
    % Draw Quad Copter
    figureHandle = uifigure("Name","Quad-Copter Simulation","Color","w","Position",[40 40 1720 780],"AutoResizeChildren","off");
    axesHandle = uiaxes(figureHandle, "Position",figureHandle.Position.*[1 1 0.6 0.9]); 
    axesHandle.Toolbar.Visible = 'off';     
    axesHandle.Interactions = [];      

    rightPanel = uipanel(figureHandle, 'Units','pixels','Position',[1200 -20 520 850],'BackgroundColor',[0.5 0.5 0.5],'BorderType','none');
    uistack(rightPanel,'top');
    %% Tabs
    tabGroup = uitabgroup(rightPanel,'Units','pixels','Position',[0 290 521 510]);   % [x y width height]

    
    tabDynamics    = uitab(tabGroup,'Title','Dynamics');
    tabEstimation   = uitab(tabGroup,'Title','Estimation');
    tabSimulation    = uitab(tabGroup,"Title",'Simulation');
    tabController   = uitab(tabGroup,'Title','Controller');
    tabWaypoints   = uitab(tabGroup,'Title','Waypoints');
    %% Dynamics tab
    uilabel(tabDynamics,'Text','Mass [kg]','Position',[20 340 80 22]);
    massField = uieditfield(tabDynamics,'numeric','Position',[110 340 80 22],'Value',USER_PARAMS.mass_kg,'Tag','MassEdit');
    
    J0 = diag([USER_PARAMS.Jx, USER_PARAMS.Jy, USER_PARAMS.Jz]);
    uilabel(tabDynamics,'Text','Inertia J [kg·m^2] (Must Be Skew Symetric!)', 'Position',[20 300 300 22]);

    JTableField = uitable(tabDynamics,'Data',J0,'ColumnEditable', true(1,3),'RowName',{'x','y','z'}, 'ColumnName',{'x','y','z'},'Position',[20 200 300 100]);

    %% Estimation tab

    rowSpacing   = 30;     % vertical space between rows
    colSpacing   = 300;    % horizontal distance between label and edit field
    startY       = 440;    % starting Y position for the first row
    labelX       = 20;     % X position of labels
    fieldWidth   = 80;
    fieldHeight  = 22;

    y = startY;
    uilabel(tabEstimation,'Text','Kalman Process Noise Estimation (deg): ', 'Position',[labelX y 260 22]);
    KalmanProcessNoiseField = uieditfield(tabEstimation,'numeric','Position',[labelX + colSpacing y fieldWidth fieldHeight], 'Value',1);
    
    y = y - rowSpacing;
    uilabel(tabEstimation,'Text','Kalman Measurement Noise Estimation (deg): ','Position',[labelX y 260 22]);
    KalmanMeasurementNoiseField = uieditfield(tabEstimation,'numeric','Position',[labelX + colSpacing y fieldWidth fieldHeight], 'Value',5);

    y = y - rowSpacing;
    uilabel(tabEstimation,'Text','Kalman Inital Angle Uncertainty (deg): ','Position',[labelX y 260 22]);
    KalmanInitalAngleUncertaintyField = uieditfield(tabEstimation,'numeric','Position',[labelX + colSpacing y fieldWidth fieldHeight], 'Value',3);
    
    y = y - rowSpacing;
    uilabel(tabEstimation,'Text','Sensor Model: Accelerometer White Noise Mag (deg): ','Position',[labelX y 300 22]);
    AccelWhiteNoiseMagField = uieditfield(tabEstimation,'numeric','Position',[labelX + colSpacing y fieldWidth fieldHeight], 'Value',3);
    
    y = y - rowSpacing;
    uilabel(tabEstimation,'Text','Sensor Model: Accelerometer Bias (deg): ', 'Position',[labelX y 300 22]);
    AccelBiasField = uieditfield(tabEstimation,'numeric', 'Position',[labelX + colSpacing y fieldWidth fieldHeight], 'Value',0);
    
    y = y - rowSpacing;
    uilabel(tabEstimation,'Text','Sensor Model: Gyroscope White Noise Mag (deg): ', 'Position',[labelX y 300 22]);
    GyroWhiteNoiseMagField = uieditfield(tabEstimation,'numeric', 'Position',[labelX + colSpacing y fieldWidth fieldHeight], 'Value',1.5);
    
    y = y - rowSpacing;
    uilabel(tabEstimation,'Text','Sensor Model: Gyroscope Bias (deg): ', 'Position',[labelX y 300 22]);
    GryoBiasField = uieditfield(tabEstimation,'numeric', 'Position',[labelX + colSpacing y fieldWidth fieldHeight], 'Value',0);

    y = y - rowSpacing;
    uilabel(tabEstimation,'Text','Enable Attitude Plotting: ', 'Position',[labelX y 200 22]);
    plotAttitudeFlagField = uiswitch(tabEstimation,'Items',{'Off','On'}, 'Value','On', 'Position',[labelX + 10  y-20  200 22]);

    y = y - 1.5*rowSpacing;
    uilabel(tabEstimation,'Text','Angle (deg) [Noise Mag Vec; Noise Frequency] (2 x N): ','Position',[labelX y 300 22]);
    NoiseMagFreqAngleField = uieditfield(tabEstimation,'text','Position',[labelX,  y-20,  fieldWidth+400,  fieldHeight],'Value','[0.01 0.02 0.03 0.04; 100 200 400 800]');
    
    y = y - 1.5*rowSpacing;
    uilabel(tabEstimation,'Text','Enable Angle Lowpass Filter: ','Position',[labelX y 260 22]);
    EnableAngleLowpassField = uiswitch(tabEstimation,'Items',{'Off','On'}, 'Value','Off', 'Position',[labelX + 10  y-20  200 22]);

    y = y - 1.5*rowSpacing;
    uilabel(tabEstimation,'Text','Angle Lowpass Cutoff (Hz): ','Position',[labelX y 260 22]);
    AngleLowpassCutoffField = uieditfield(tabEstimation,'numeric','Position',[labelX + colSpacing y fieldWidth fieldHeight], 'Value',15);

    %% Simulation tab
    ySim = 440;
    colSpacing = 200;
    uilabel(tabSimulation,'Text','Simulation Hz: ','Position',[labelX ySim 100 22]);
    simHzField = uieditfield(tabSimulation,'numeric','Position',[labelX + colSpacing ySim 80 22], 'Value',150,'Tag','Hz');
    
    ySim = ySim - rowSpacing;
    uilabel(tabSimulation,'Text','Display Hz: ','Position',[labelX ySim 80 22]);
    displayHzField = uieditfield(tabSimulation,'numeric','Position',[labelX + colSpacing ySim fieldWidth fieldHeight],'Value',30,'Tag','display Hz');
    
    ySim = ySim - rowSpacing;
    uilabel(tabSimulation,'Text','Quad Model Scale: ','Position',[labelX ySim 260 22]);
    QuadScaleField = uieditfield(tabSimulation,'numeric','Position',[labelX + colSpacing ySim fieldWidth fieldHeight], 'Value',USER_PARAMS.quadScale);

    ySim = ySim - rowSpacing;
    uilabel(tabSimulation,'Text','Prop Spin Speed (rpm): ','Position',[labelX ySim 260 22]);
    PropSpinRPMField = uieditfield(tabSimulation,'numeric','Position',[labelX + colSpacing ySim fieldWidth fieldHeight], 'Value',USER_PARAMS.propSpinRPM);

    %% Controller tab

    y = startY; 
    y = y - rowSpacing;
    uilabel(tabController,'Text','Outer Loop: Angular Velocity (deg/s): ', 'Position',[labelX y 300 22]);
    OuterLoopAngularVelocityField = uieditfield(tabController,'numeric', 'Position',[labelX + colSpacing y fieldWidth fieldHeight], 'Value',90);
    
    y = y - rowSpacing;
    uilabel(tabController,'Text','Outer Loop: Acceleration (mps2): ', 'Position',[labelX y 300 22]);
    OuterLoopAccelerationField = uieditfield(tabController,'numeric', 'Position',[labelX + colSpacing y fieldWidth fieldHeight], 'Value',100);

    y = y - rowSpacing;
    uilabel(tabController,'Text','Angle Cmd Smoothing Filter (sec): ','Position',[labelX y 260 22]);
    AngleCommandFilterTauField = uieditfield(tabController,'numeric','Position',[labelX + colSpacing y fieldWidth fieldHeight],'Value',USER_PARAMS.AngleCommandFilterTau);

    y = y - 1.5*rowSpacing;
    uilabel(tabController,'Text','Enable Angle Controller (False: Velocity Controller): ', 'Position',[labelX y 300 22]);
    EnableAngleControlField = uiswitch(tabController,'Items',{'Off','On'}, 'Value','On', 'Position',[labelX + 10  y-20  200 22]);

    y = y - 1.5*rowSpacing;
    uilabel(tabController,'Text','Enable Adaptive Controller (False: Fixed Gain): ', 'Position',[labelX y 300 22]);
    EnableAdaptiveControlField = uiswitch(tabController,'Items',{'Off','On'}, 'Value','On', 'Position',[labelX + 10  y-20  200 22]);

    % DMAC settings
    colSpacing = colSpacing*1.2;
    y = y - 2*rowSpacing;
    uilabel(tabController,'Text',' DMAC Tunnable Optimization Parameters: ','Position',[labelX y 260 22],'FontWeight','bold');

    y = y - rowSpacing;
    uilabel(tabController,'Text','DMAC: Q diag (9 elements):','Position',[labelX y 260 22]);
    DMAC_Qdiag_Field = uieditfield(tabController,'text','Position',[labelX + colSpacing, y, fieldWidth+150, fieldHeight],'Value','[1e2 1e2 1e2 1e0 1e0 1e0 1e1 1e1 1e1]');
    
    y = y - rowSpacing;
    uilabel(tabController,'Text','DMAC: R diag (3 elements):','Position',[labelX y 260 22]);
    DMAC_Rdiag_Field = uieditfield(tabController,'text','Position',[labelX + colSpacing, y, fieldWidth+20, fieldHeight],'Value','[1e-3 1e-3 1e-3]');
    
    y = y - rowSpacing;
    uilabel(tabController,'Text','DMAC: R0: ','Position',[labelX y 260 22]);
    DMAC_R0_Field = uieditfield(tabController,'numeric','Position',[labelX + colSpacing, y, fieldWidth, fieldHeight],'Value',1e-5);
    
    y = y - rowSpacing;
    uilabel(tabController,'Text','DMAC: Forgetting Factor half-life (sec): ','Position',[labelX y 260 22]);
    DMAC_lambda_Field = uieditfield(tabController,'numeric','Position',[labelX + colSpacing, y, fieldWidth, fieldHeight],'Value',1.5);
    
    %% Waypoint Tab
    y = startY;

    % Enable waypoint navigation switch
    uilabel(tabWaypoints,'Text','Enable Waypoint Navigation: ','Position',[labelX y 260 22]);
    EnableWaypointField = uiswitch(tabWaypoints,'Items',{'Off','On'},'Value','Off','Position',[labelX + 10, y-20, 200, 22]);

    y = y - 2*rowSpacing;
    uilabel(tabWaypoints,'Text','Waypoints [3 x N], [East List; North List; Up List]:','Position',[labelX y 400 22]);

    WaypointDataField = uieditfield(tabWaypoints,'text','Position',[labelX,  y-20,  fieldWidth+300,  fieldHeight],'Value','[10 20 20 10; 10 10 20 20; 10 20 30 40]');

    y = y - 2*rowSpacing;
    uilabel(tabWaypoints,'Text','Waypoint Cruise Speed (m/s):','Position',[labelX y 260 22]);
    WaypointCruiseSpeedField = uieditfield(tabWaypoints,'numeric','Position',[labelX + colSpacing, y, fieldWidth, fieldHeight],'Value',USER_PARAMS.waypointCruiseSpeed);

    y = y - rowSpacing;
    uilabel(tabWaypoints,'Text','Waypoint Slow Radius (m):','Position',[labelX y 260 22]);
    WaypointSlowRadiusField = uieditfield(tabWaypoints,'numeric', 'Position',[labelX + colSpacing, y, fieldWidth, fieldHeight],'Value',USER_PARAMS.waypointSlowRadius);
    %% pack all the stuff into the objects struct:
    % dynamics
    OBJs.GUI_Params.massField = massField;
    OBJs.GUI_Params.JTableField = JTableField;
    % estimation 
    OBJs.GUI_Params.KalmanProcessNoiseField = KalmanProcessNoiseField;
    OBJs.GUI_Params.KalmanMeasurementNoiseField = KalmanMeasurementNoiseField;
    OBJs.GUI_Params.KalmanInitalAngleErrorField = KalmanInitalAngleUncertaintyField;
    OBJs.GUI_Params.PlotAttitudeField = plotAttitudeFlagField;
    OBJs.GUI_Params.AccelWhiteNoiseMagField = AccelWhiteNoiseMagField;
    OBJs.GUI_Params.GyroWhiteNoiseMagField = GyroWhiteNoiseMagField;
    OBJs.GUI_Params.AccelBiasField = AccelBiasField;
    OBJs.GUI_Params.GryoBiasField = GryoBiasField;
    OBJs.GUI_Params.NoiseMagFreqAngleField = NoiseMagFreqAngleField;
    OBJs.GUI_Params.NoiseMagFreqAngleField = NoiseMagFreqAngleField;
    OBJs.GUI_Params.EnableAngleLowpassField = EnableAngleLowpassField;
    OBJs.GUI_Params.AngleLowpassCutoffField = AngleLowpassCutoffField;
    % simulation
    OBJs.GUI_Params.simHzField = simHzField;
    OBJs.GUI_Params.displayHzField = displayHzField;
    OBJs.GUI_Params.QuadScaleField = QuadScaleField;
    OBJs.GUI_Params.PropSpinRPMField = PropSpinRPMField;
    % control
    OBJs.GUI_Params.OuterLoopAngularVelocityField = OuterLoopAngularVelocityField;
    OBJs.GUI_Params.OuterLoopAccelerationField = OuterLoopAccelerationField;
    OBJs.GUI_Params.AngleCommandFilterTauField = AngleCommandFilterTauField;
    OBJs.GUI_Params.EnableAdaptiveControlField = EnableAdaptiveControlField;
    OBJs.GUI_Params.EnableAngleControlField = EnableAngleControlField;
    OBJs.GUI_Params.DMAC_Qdiag_Field = DMAC_Qdiag_Field;
    OBJs.GUI_Params.DMAC_Rdiag_Field = DMAC_Rdiag_Field;
    OBJs.GUI_Params.DMAC_R0_Field = DMAC_R0_Field;
    OBJs.GUI_Params.DMAC_lambda_Field = DMAC_lambda_Field;
    % waypoint
    OBJs.GUI_Params.EnableWaypointField = EnableWaypointField;
    OBJs.GUI_Params.WaypointDataField  = WaypointDataField;
    OBJs.GUI_Params.WaypointCruiseSpeedField = WaypointCruiseSpeedField;
    OBJs.GUI_Params.WaypointSlowRadiusField = WaypointSlowRadiusField;
    %% Start, Pause, Stop Buttons
   
    uibutton(rightPanel,'Text','Start','Position',[40 40 140 40], 'ButtonPushedFcn',@(btn,~) startSimCB(figureHandle));
    
    uibutton(rightPanel,'Text','Pause','Position',[200 40 140 40],'ButtonPushedFcn',@(btn,~) stopSimCB(figureHandle));
    
    uibutton(rightPanel,'Text','Restart','Position',[360 40 140 40],'ButtonPushedFcn',@(btn,~) restartSimCB(figureHandle));

    %% Draw
    grid(axesHandle,'on'); box(axesHandle,'on'); axis(axesHandle,'equal');
    xlim(axesHandle, [-USER_PARAMS.xLim,USER_PARAMS.xLim]); ylim(axesHandle, [-USER_PARAMS.xLim,USER_PARAMS.xLim]); zlim(axesHandle, [-0.1,250]);
    zPlane = [-0.01 -0.01 -0.01 -0.01]; 
    xPlane = [USER_PARAMS.xLim, -USER_PARAMS.xLim, -USER_PARAMS.xLim, USER_PARAMS.xLim];
    yPlane = [USER_PARAMS.xLim, USER_PARAMS.xLim, -USER_PARAMS.xLim, -USER_PARAMS.xLim];
    pH = patch(axesHandle, "XData", xPlane, "YData",yPlane, "ZData",zPlane);
    pH.FaceColor = [0.7 0.7 0.7];
    
    camup(axesHandle, [0 0 1]);
    axesHandle.Projection = "perspective";    

    scale = USER_PARAMS.quadScale;
    %  Quad geometry
    armLength = 0.35*scale;      
    armWidth = 0.04*scale;      
    armHeight = 0.02*scale;      
    motWidth = 0.06*scale;      
    motHeight = 0.03*scale;     
    propRadius = 0.13*scale;    
        
    % Arm boxes in body frame centered at origin, lying along X and Y
    [VarmX,Fbox] = boxVerts([-armLength/2, -armWidth/2, -armHeight/2], [armLength/2, armWidth/2, armHeight/2]); % X arm
    [VarmY,~]    = boxVerts([-armWidth/2, -armLength/2, -armHeight/2], [armWidth/2, armLength/2, armHeight/2]); % Y arm
    
    % Motor centers (X+, X-, Y+, Y-)
    motors = [ armLength/2, 0, 0;
              -armLength/2, 0, 0;
               0,  armLength/2, 0;
               0, -armLength/2, 0];
    
    % Transforms 
    tQuad = hgtransform('Parent',axesHandle);                    % body pose
    % Arms (as single patches under body transform)
    hArmX = patch('Faces',Fbox,'Vertices',VarmX,'FaceColor',[0.1 0.1 0.1],'EdgeColor','none','Parent',tQuad);
    hArmY = patch('Faces',Fbox,'Vertices',VarmY,'FaceColor',[0.2 0 0],'EdgeColor','none','Parent',tQuad);

    tProp = gobjects(4,1);
    hMot  = gobjects(4,1);
    hBlade = gobjects(4,4);    

    for i=1:4
        tProp(i) = hgtransform('Parent', tQuad);  % local node at motor i
    
        % motor as a cylinder centered at motors(i,:)
        motR = 0.5*motWidth;              
        nSides = 40;
        [Xc,Yc,Zc] = cylinder(motR, nSides);   
        Zc = Zc * motHeight;                        
        hMot(i) = surf(Xc + motors(i,1), Yc + motors(i,2), Zc + motors(i,3), 'Parent', tQuad, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor','none', 'SpecularStrength',0.3, 'DiffuseStrength',0.8, 'AmbientStrength',0.3);
        Rb = propRadius;  Wb = 0.05*scale;  zOff = motHeight*1.01; 
        % Blade 1 along +X
        B1 = [ 0   -Wb/2  zOff;
               Rb -Wb/2  zOff;
               Rb  Wb/2  zOff;
               0   Wb/2  zOff ];
    
        B2 = (R_from_phi_theta_psi(0,0,90)*B1.').';
        B3 = (R_from_phi_theta_psi(0,0,180)*B1.').';
        B4 = (R_from_phi_theta_psi(0,0,270)*B1.').';
    
        % Create blades as two quads under the prop transform
        hBlade(i,1) = patch('Vertices',B1,'Faces',[1 2 3 4],'FaceColor',[1 0 0.7], 'EdgeColor','none', 'Parent', tProp(i));
        hBlade(i,2) = patch('Vertices',B2,'Faces',[1 2 3 4],'FaceColor',[1 0 0.7], 'EdgeColor','none', 'Parent', tProp(i));
        hBlade(i,3) = patch('Vertices',B3,'Faces',[1 2 3 4],'FaceColor',[1 0 0.7], 'EdgeColor','none', 'Parent', tProp(i));
        hBlade(i,4) = patch('Vertices',B4,'Faces',[1 2 3 4],'FaceColor',[1 0 0.7], 'EdgeColor','none', 'Parent', tProp(i));
        % Place prop transform at motor center in body frame
        set(tProp(i),'Matrix', makehgtform('translate', motors(i,:)));
    end
    
        
    % Shadow setup
    OBJs.DrawShadow = USER_PARAMS.DrawShadow;
    % sun position
    OBJs.LightPosWorld = [200, -300, 1000];  
    % collect body-frame vertices for shadow hull
    V_body = [];
    V_body = [V_body; VarmX; VarmY];                 
    % motor cylinder vertices in body frame
    for i = 1:4
        motR = 0.5*motWidth; nSides = 40;
        [Xc,Yc,Zc] = cylinder(motR, nSides); 
        Zc = Zc * motHeight;
        Vc = [Xc(:)+motors(i,1), Yc(:)+motors(i,2), Zc(:)+motors(i,3)];
        V_body = [V_body; Vc];
    end

    V_body = [V_body; B1; B2; B3; B4];
    OBJs.V_body = V_body; 
    
    % shadow patch on ground
    hShadow = patch(axesHandle,'XData', nan, 'YData', nan, 'ZData', nan,'FaceColor', [0 0 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

    OBJs.hShadow = hShadow;
    OBJs.tProp = tProp;
    OBJs.motors = motors;
    OBJs.tQuad = tQuad;

    %% attitude plot
    positionMain = axesHandle.Position;        % [x y w h] in pixels

    insetWidth = 0.25*positionMain(3);             % 25% of main width
    insetHieght = 0.25*positionMain(4);             % 25% of main height

    insetPos = [positionMain(1) - 10, positionMain(2) - 10, insetWidth, insetHieght];

    attitudeAxes = uiaxes(figureHandle,'Position', insetPos,'Box','on','Color','w');
    xlabel(attitudeAxes,'k');
    ylabel(attitudeAxes,'angle (deg)');
    hold(attitudeAxes,'on');

    attitudeAxes.YLimMode = 'auto';
    attitudeAxes.Box = 'on';
    attitudeAxes.FontSize = 8;
    title(attitudeAxes,'Attitude [red: pitch, blue: roll, - truth, -- estimation');
    OBJs.attitudeAxes = attitudeAxes;

    %% position plot
    insetPositionPositionPlot = [positionMain(1) + positionMain(3) - insetWidth + 100, positionMain(2) - 10, insetWidth, insetHieght];

    positionAxes = uiaxes(figureHandle,'Position', insetPositionPositionPlot,'Box','on','Color','w');
    xlabel(positionAxes,'east (m)');
    ylabel(positionAxes,'north (m)');
    hold(positionAxes,'on');

    positionAxes.YLimMode = 'auto';
    positionAxes.Box = 'on';
    positionAxes.FontSize = 8;
    title(positionAxes,'Top Down Position (east vs north)');
    OBJs.positionAxes = positionAxes;
    %% FFT Plot

    insetPositionFFT = [positionMain(1) - 10, positionMain(2) + positionMain(4) - insetHieght + 10, insetWidth, insetHieght];

    fftAxes = uiaxes(figureHandle,'Position', insetPositionFFT,'Box','on','Color','w');
    xlabel(fftAxes,'f (Hz)');
    ylabel(fftAxes,'|Theta|');
    fftAxes.YLimMode = 'auto';
    fftAxes.Box = 'on';
    fftAxes.FontSize = 8;
    title(fftAxes,'FFT of \theta_{est}');
    OBJs.fftAxes = fftAxes;

end

function OBJs = updateQuadPlot(iter, axesHandle, states, dt, OBJs)
    persistent spin
    if iter == 1
        spin  = 0;
    end

    tQuad = OBJs.tQuad;
    tProp = OBJs.tProp;
    motors = OBJs.motors;

    spinGainRPM = OBJs.GUI_Params.PropSpinRPMField.Value;
    spinGain = spinGainRPM*0.10472;
    spinDir = [-1 1 1 -1];              

    phi = states(1); theta = states(2); psi = -states(3);
    x = -states(8); y = -states(7); z = states(9);
   
    R = R_from_phi_theta_psi(phi, theta, psi);
    T = makehgtform('translate', [x y z]) * [R [0;0;0]; 0 0 0 1];
    set(tQuad,'Matrix',T);

    if isfield(OBJs,'DrawShadow') && OBJs.DrawShadow
        updateShadow(OBJs.hShadow, OBJs.V_body, OBJs.LightPosWorld, x, y, z, R);
    end

    spin = spin + spinGain*dt;


    for i=1:4
        set(tProp(i),'Matrix', makehgtform('translate', motors(i,:)) * makehgtform('zrotate', spinDir(i)*spin));
    end

    persistent camPos camTgt PlotAttitudeBoolPrev OriginalPosition
    if iter == 1
        camPos = [x-10, y-10, z+10];
        camTgt = [x, y, z];
        axis(axesHandle, 'vis3d');
        daspect(axesHandle, [1 1 1]);
        camproj(axesHandle, 'perspective');
        camva(axesHandle, 60);
        camup(axesHandle, [0 0 1]);
        PlotAttitudeBoolPrev = 'on';
        OriginalPosition = [];
    end
    R = R_from_phi_theta_psi(0, 0, psi);
    forward = R * [1;0;0];
    offset = 15 * forward + [0; 0; 10];
    desiredCam = [x; y; z] + offset;
    alpha = exp(-dt/0.015);  % smoothing factor

    camPos = alpha * camPos + (1-alpha) * desiredCam';
    camTgt = alpha * camTgt + (1-alpha) * [x, y, z];

    campos(axesHandle, camPos);
    camtarget(axesHandle, camTgt);   

    
    if isfield(OBJs,'attitudeAxes') && isvalid(OBJs.attitudeAxes)
        attitudeAxes = OBJs.attitudeAxes;
        if isfield(OBJs,'AttitudeLog')
            AttitudeLog = OBJs.AttitudeLog;
        else
            AttitudeLog = [];
        end

        if ~isempty(AttitudeLog)
            maxN = 300;
            n = size(AttitudeLog,2);
            idx = max(1,n-maxN):n;
            Aplot = AttitudeLog(:,idx);

            cla(attitudeAxes)
            plot(attitudeAxes, idx, Aplot(1,:)*180/pi, 'r',  'LineWidth',1); % phi true
            plot(attitudeAxes, idx, Aplot(2,:)*180/pi, 'r--','LineWidth',1);% phi est
            plot(attitudeAxes, idx, Aplot(3,:)*180/pi, 'b',  'LineWidth',1); % theta true
            plot(attitudeAxes, idx, Aplot(4,:)*180/pi, 'b--','LineWidth',1); % theta est
            
            if length(idx) > 1
                attitudeAxes.XLim = [idx(1) idx(end)];
            end
        end

         if isempty(OriginalPosition)
            OriginalPosition = attitudeAxes.Position;
         end

        % only do it if its new to save cpu 
        if ~strcmp(PlotAttitudeBoolPrev,OBJs.GUI_Params.PlotAttitudeField.Value)
            if strcmp(OBJs.GUI_Params.PlotAttitudeField.Value,'On')
                attitudeAxes.Visible = 'on';
                % disp('switched to on')
           
                for I = 1:length(attitudeAxes.Children)
                    attitudeAxes.Children(I).Visible = 'on';
                    attitudeAxes.Position = OriginalPosition;
                end
            else
                attitudeAxes.Visible = 'off';
                attitudeAxes.Position = [0 0 0 0];
            end
        end
        PlotAttitudeBoolPrev = OBJs.GUI_Params.PlotAttitudeField.Value;
    end

    if isfield(OBJs,'positionAxes') && isvalid(OBJs.positionAxes)
        positionAxes = OBJs.positionAxes;
        if isfield(OBJs,'PositionLog')
            PositionLog = OBJs.PositionLog;
        else
            PositionLog = [];
        end

        if ~isempty(PositionLog)
            maxN = 1000;
            n = size(PositionLog,2);
            idx = max(1,n-maxN):n;
            
            east  = -PositionLog(1,idx);  
            north =  PositionLog(2,idx);

            cla(positionAxes)
            hold(positionAxes,'on');

            % defaults in case waypoint GUI fields are missing
            EnableWaypointFlag = 'Off';
            wp_east  = [];
            wp_north = [];

            % draw waypoints circles if enabled
            if isfield(OBJs,'GUI_Params') && isfield(OBJs.GUI_Params,'EnableWaypointField')
                EnableWaypointFlag = OBJs.GUI_Params.EnableWaypointField.Value;
                if strcmp(EnableWaypointFlag,'On')
                    wpString = OBJs.GUI_Params.WaypointDataField.Value;
                    waypointList = str2num(wpString);

                    if ~isempty(waypointList) && size(waypointList,1) == 3
                        wp_east  = -waypointList(1,:);  
                        wp_north =  waypointList(2,:);
                        plot(positionAxes, wp_east, wp_north, 'ko','MarkerSize',6, 'LineWidth',1); 
                    end
                end
            end

            % draw trajectory
            plot(positionAxes, east, north, 'r', 'LineWidth',2);

            if length(idx) > 1
                allX = [east,  wp_east];
                allY = [north, wp_north];

                margin = 5;
                positionAxes.XLim = [min(allX)-margin, max(allX)+margin];
                positionAxes.YLim = [min(allY)-margin, max(allY)+margin];
            end
        end
    end




     if isfield(OBJs,'fftAxes') && isvalid(OBJs.fftAxes) && isfield(OBJs,'fAxis') && isfield(OBJs,'P1Theta')
    
        fftAxes = OBJs.fftAxes;
        cla(OBJs.fftAxes);

        plot(fftAxes, OBJs.fAxis, OBJs.P1Theta, 'k-','LineWidth',1);
        xlim(fftAxes,[0 max(OBJs.fAxis)]);
        fftAxes.XLimMode = 'manual'; 
        fftAxes.YLimMode = 'auto';
    end
end

%% Local Helpers
function velocityCommand = getWaypointVelocityCommandFromWaypointList(resetFlag, waypointList, stateVector, USER_PARAMS)
    persistent wayPointIndex

    if isempty(wayPointIndex) || resetFlag
        wayPointIndex = 1;
    end
    position = stateVector(7:9,1); % (ENU) (3x1);
    wayPointPosition = waypointList(:,wayPointIndex);

    positionError = wayPointPosition - position;

    cruiseSpeed = USER_PARAMS.waypointCruiseSpeed; %proportional gain
    
    % position error is normed so the velocity components stay bounded
    R = norm(positionError);
    velocityCommandCruise = cruiseSpeed * positionError/R; 
    
    % start slowing when close to waypoint
    slowRadius = USER_PARAMS.waypointSlowRadius;
    if R < slowRadius
        velocityCommand = velocityCommandCruise * R/slowRadius;
    else
        velocityCommand = velocityCommandCruise;
    end

    if R < 1 %4 is a magic number i made up, switched if distance is 4x smaller than the slow radius
        wayPointIndex = wayPointIndex + 1;
        if length(waypointList(1,:)) < wayPointIndex
            wayPointIndex = 1;
        end
    end
end

function moments = innerLoopController(iter, dt, states, angle_command, params)
    persistent K discreteErrorIntegral
    if iter == 1 || isempty(discreteErrorIntegral)
        % A_cont = zeros(6,6); 
        % A_cont(1:3,4:6) = eye(3);
        % B_cont = zeros(6,3);
        % B_cont(4:6,:) = params.J;
        % 
        % C = zeros(3,6);
        % C(:,1:3) = eye(3);
        % 
        % D = [];
        % 
        % system = ss(A_cont,B_cont,C,D);
        % 
        % system_discrete = c2d(system,dt);
        % [Aa, Ba] = augmentMatriciesDiscrete(system_discrete.A, system_discrete.B, system_discrete.C);
        % 
        % desiredPoles = exp(linspace(-20,-10,9).*dt);
        % K = -place(Aa,Ba, desiredPoles);

        discreteErrorIntegral = zeros(3,1);
        % decided to hard code for simplicity
        K = 1e3.*[  -0.981937335066258   0.008465455222118   0.003067159723532  -0.037323909691861   0.000656492127151  -0.000209610752826   0.073154504878000   0.000578296421453   0.000097162124784;
                     0.032534297939173  -1.030398816832853   0.013708568293784   0.000406953635326  -0.041006539208563   0.000390729276341  -0.000911357929482   0.073241772607058  -0.000309567882199;
                     0.004725624329662   0.005453323431547  -0.745105950765370   0.000578002938297   0.000294179965578  -0.026545953391705  -0.000124764184386   0.000129940083088   0.059677716439461];
    end
    angles = states(1:3,1);
    innerLoopStates = states(1:6,1);
    discreteErrorIntegral = discreteErrorIntegral + (angle_command-angles);
    moments = K*[innerLoopStates; discreteErrorIntegral];   
end

function [Aa, Ba] = augmentMatriciesDiscrete(A, B, C)

    [ny, ~] = size(C);
    [~, nu] = size(B);
    nx = size(A,1);

    D = zeros(ny,nu);
 
    Aa = [A, zeros(nx,ny);
         -C,   eye(ny,ny)];
    Ba = [B; -D];
end

function xk1 = RK4(func, dT, x, moments, force, params)
  k1 = func(x,moments, force, params);
  k2 = func(x+dT/2*k1,moments, force, params);
  k3 = func(x+dT/2*k2,moments, force, params);
  k4 = func(x+dT*k3,  moments, force, params);
  xk1 = x + dT/6 * (k1 + 2*k2 + 2*k3 + k4); 
end

function DMAC_STATES = initDMAC_Settings(dt)
    DMAC_STATES.UsePlaceAndLQR = false;
    DMAC_STATES.Q = diag([1e2 1e2 1e2 1e0 1e0 1e0 1e1 1e1 1e1]);
    DMAC_STATES.R = diag([1e-3,1e-3,1e-3]);
    DMAC_STATES.R0 = 1e-5;
    DMAC_STATES.C = [eye(3), zeros(3)];
       
    halfLife = 1.5;
    DMAC_STATES.lambda = 0.5^(dt/halfLife);
    DMAC_STATES.acheivedPoles = nan(9,1);
    DMAC_STATES.K = nan;
end

function [phi_hat, theta_hat] = KalmanFilterRollPitch(dt, gyro_pqr, accel_xyz, reset, initAngleUncertainty, measurementNoise, processNoise)

    persistent phi_est theta_est P_phi P_theta

    if isempty(phi_est) || reset
        % Use accel to get initial attitude
        ax = accel_xyz(1);
        ay = accel_xyz(2);
        az = accel_xyz(3);
        denom = max(sqrt(ay*ay + az*az), 1e-9);

        phi_meas   = atan2(-ax, denom);
        theta_meas = atan2( ay, az );

        phi_est   = phi_meas;
        theta_est = theta_meas;

        % Initial angle uncertainty
        sigma0 = deg2rad(initAngleUncertainty);
        P_phi   = sigma0^2;
        P_theta = sigma0^2;
    end

    p = gyro_pqr(1);
    q = gyro_pqr(2);
    r = gyro_pqr(3);

    ax = accel_xyz(1);
    ay = accel_xyz(2);
    az = accel_xyz(3);
   
    % [phi_dot; theta_dot; psi_dot] = S_Phi_Theta_inv(phi,theta) * [p;q;r]
    S = S_Phi_Theta_inv(phi_est, theta_est);
    eul_dot = S * [p; q; r];

    phi_dot   = eul_dot(1);
    theta_dot = eul_dot(2);

    % Predict angles
    phi_pred   = phi_est   + dt * phi_dot;
    theta_pred = theta_est + dt * theta_dot;

    % processNoise = deg2rad(1);              % rad/s
    processNoise = deg2rad(processNoise);
    Q_phi   = (processNoise)^2;
    Q_theta = (processNoise)^2;

    P_phi   = P_phi   + (dt^2) * Q_phi;
    P_theta = P_theta + (dt^2) * Q_theta;

    phi_meas = atan2(ay, sqrt(ax^2 + az^2));
    theta_meas = -atan2(ax, sqrt(ay^2 +az^2));

    % Measurement noise
    % r_meas = deg2rad(5);
    r_meas = deg2rad(measurementNoise);
    R_phi   = r_meas^2;
    R_theta = r_meas^2;

    S_phi = P_phi + R_phi;
    K_phi = P_phi / S_phi;

    phi_est = phi_pred + K_phi * (phi_meas - phi_pred);
    P_phi   = (1 - K_phi) * P_phi;

    S_theta = P_theta + R_theta;
    K_theta = P_theta / S_theta;

    theta_est = theta_pred + K_theta * (theta_meas - theta_pred);
    P_theta   = (1 - K_theta) * P_theta;

    phi_est   = wrapToPi(phi_est);
    theta_est = wrapToPi(theta_est);

    phi_hat   = phi_est;
    theta_hat = theta_est;
    assignin('base','phi_hat_last',phi_hat);
    assignin('base','theta_hat_last',theta_hat);
end

function [filtered_phi_hat, filtered_theta_hat] = lowpassFilter(reset, phi_hat, theta_hat, cutOffFreqHz, dt)
    persistent internalFilteringStatePhi internalFilteringStateTheta A B C D

    if reset || isempty(internalFilteringStatePhi)
        s = tf('s');
        cutOffFreqRPS = 2*pi*cutOffFreqHz;
        TF_LP = 1/(1+s/cutOffFreqRPS);
        TF_LP_Discrete = c2d(TF_LP, dt);
        [A, B, C, D] = tf2ss(TF_LP_Discrete.Numerator{1}, TF_LP_Discrete.Denominator{1});
        internalFilteringStatePhi = phi_hat;
        internalFilteringStateTheta = theta_hat;
    end

    internalFilteringStatePhi = A*internalFilteringStatePhi + B*phi_hat;
    filtered_phi_hat = C*internalFilteringStatePhi + D*phi_hat;

    internalFilteringStateTheta = A*internalFilteringStateTheta + B*theta_hat;
    filtered_theta_hat = C*internalFilteringStateTheta + D*theta_hat;

    % testing script to test on step function
    % x = [zeros(75,1); ones(75,1)];
    % x = x + 0.01*randn(150,1);
    % for I = 1:150
    %     if I == 1
    %         y(I) = lowpassFilter(true, x(I), x(I), 10, 1/150);
    %     else
    %         y(I) = lowpassFilter(false, x(I), x(I), 10, 1/150);
    %     end
    % end
    % figure; hold on;
    % plot(1/150*(1:150),x);
    % plot(1/150*(1:150),y);
end

function S = S_Phi_Theta_inv(Phi,Theta)
    S = [1 sin(Phi)*tan(Theta) cos(Phi)*tan(Theta);
        0 cos(Phi) -sin(Phi);
        0 sin(Phi)*sec(Theta) cos(Phi)*sec(Theta)];
end
%% Plotting helpers
function R = R_from_phi_theta_psi(phi,theta,psi)
    % phi = pitch about Y, theta = roll about X, psi = yaw about Z
    Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    Ry = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
    Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    R  = Rz*Ry*Rx;
end

function [V,F] = boxVerts(lo,hi)
    % Generated via chat GPT for making the shadow
    % axis-aligned box from lo=[x y z] to hi=[x y z]
    x = [lo(1) hi(1)]; y = [lo(2) hi(2)]; z = [lo(3) hi(3)];
    V = [x(1) y(1) z(1);
         x(2) y(1) z(1);
         x(2) y(2) z(1);
         x(1) y(2) z(1);
         x(1) y(1) z(2);
         x(2) y(1) z(2);
         x(2) y(2) z(2);
         x(1) y(2) z(2)];
    F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
end

function updateShadow(hShadow, V_body, lightPosWorld, x, y, z, R)
    % transform body verts to world
    % generated by chat GPT
    Vw = (R * V_body.').';
    Vw = Vw + repmat([x y z], size(Vw,1), 1);
    % project to z=0
    Vsh = projectShadowToZ0(Vw, lightPosWorld);
    try
        K = convhull(Vsh(:,1), Vsh(:,2));
        Xs = Vsh(K,1); Ys = Vsh(K,2); Zs = zeros(size(K));
    catch
        Xs = Vsh(:,1); Ys = Vsh(:,2); Zs = zeros(size(Vsh,1),1);
    end
    set(hShadow, 'XData', Xs, 'YData', Ys, 'ZData', Zs);
end

function Vshadow = projectShadowToZ0(Vworld, Lpos) 
    % generated by chat GPT
    Lx = Lpos(1); Ly = Lpos(2); Lz = Lpos(3); X = Vworld(:,1); 
    Y = Vworld(:,2); Z = Vworld(:,3); 
    denom = (Lz - Z); 
    denom(abs(denom) < 1e-9) = sign(denom(abs(denom) < 1e-9))*1e-9; 
    Sx = Lz .* (X - Lx) ./ denom + Lx; 
    Sy = Lz .* (Y - Ly) ./ denom + Ly; 
    Vshadow = [Sx, Sy, zeros(size(Sx))]; 
end
%% Keyboard input detection helpers and call back helpers
function [userCommands, velocityCommands, psi_cmd] = getKeysAndVelocityCommands(userCommands, stateVector, USER_PARAMS, dt)
        % update commands based on user's keyboard
        % Returns true/false if key is down

        % forward and back
        if evalin('base','exist(''KEYS_W'',''var'') && KEYS_W')
            userCommands.velocityY = userCommands.velocityY + USER_PARAMS.acceleration*dt;
        elseif evalin('base','exist(''KEYS_S'',''var'') && KEYS_S')
            userCommands.velocityY = userCommands.velocityY - USER_PARAMS.acceleration*dt;
        else
            userCommands.velocityY = 0;
        end
        
        % left and right
        if evalin('base','exist(''KEYS_D'',''var'') && KEYS_D')
            userCommands.velocityX = userCommands.velocityX - USER_PARAMS.acceleration*dt;
        elseif evalin('base','exist(''KEYS_A'',''var'') && KEYS_A')
            userCommands.velocityX = userCommands.velocityX + USER_PARAMS.acceleration*dt;
        else
            userCommands.velocityX = 0;
        end
               
        % rotation
        if evalin('base','exist(''KEYS_E'',''var'') && KEYS_E')
            userCommands.psi = userCommands.psi + USER_PARAMS.angular_velocity*dt;
        end
        if evalin('base','exist(''KEYS_Q'',''var'') && KEYS_Q')
            userCommands.psi = userCommands.psi - USER_PARAMS.angular_velocity*dt;
        end
        
        % up and down
        if evalin('base','exist(''KEYS_SHIFT'',''var'') && KEYS_SHIFT')
            userCommands.velocityZ = userCommands.velocityZ + USER_PARAMS.acceleration*dt;
        elseif evalin('base','exist(''KEYS_CONTROL'',''var'') && KEYS_CONTROL')
            userCommands.velocityZ = userCommands.velocityZ - USER_PARAMS.acceleration*dt;
        else
            userCommands.velocityZ = 0;
        end
        
        % safety clamping
        userCommands.velocityX = clip(userCommands.velocityX, -USER_PARAMS.maxVelocity, USER_PARAMS.maxVelocity);
        userCommands.velocityY = clip(userCommands.velocityY, -USER_PARAMS.maxVelocity, USER_PARAMS.maxVelocity);
        userCommands.velocityZ = clip(userCommands.velocityZ, -USER_PARAMS.maxVelocity, USER_PARAMS.maxVelocity);
    
        % for now, preventing psi to equal pi
        userCommands.psi = clip(wrapToPi(userCommands.psi),-0.995*pi,0.995*pi);
        psi_cmd = userCommands.psi;
        
        % user commands are in body fixed frame
        phi = stateVector(1); theta = stateVector(2); psi = stateVector(3);
        R_body_to_world = R_from_phi_theta_psi(phi, theta, psi);
        
        vel_body = [userCommands.velocityX; userCommands.velocityY; userCommands.velocityZ];
        vel_world = R_body_to_world * vel_body;
        
        velocityCommands = vel_world;
end

function [userCommands, angleCommands, force] = getKeysAndAngleCommands(userCommands, stateVector, USER_PARAMS, dt)
      % forward and back
        if evalin('base','exist(''KEYS_W'',''var'') && KEYS_W')
            userCommands.pitch = userCommands.pitch - USER_PARAMS.acceleration*dt;
        elseif evalin('base','exist(''KEYS_S'',''var'') && KEYS_S')
            userCommands.pitch = userCommands.pitch + USER_PARAMS.acceleration*dt;
        else
            userCommands.pitch = 0;
        end
        
        % left and right
        if evalin('base','exist(''KEYS_D'',''var'') && KEYS_D')
            userCommands.roll = userCommands.roll - USER_PARAMS.acceleration*dt;
        elseif evalin('base','exist(''KEYS_A'',''var'') && KEYS_A')
            userCommands.roll = userCommands.roll + USER_PARAMS.acceleration*dt;
        else
            userCommands.roll = 0;
        end
               
        % rotation
        if evalin('base','exist(''KEYS_Q'',''var'') && KEYS_Q')
            userCommands.psi = userCommands.psi - USER_PARAMS.angular_velocity*dt;
        end
        if evalin('base','exist(''KEYS_E'',''var'') && KEYS_E')
            userCommands.psi = userCommands.psi + USER_PARAMS.angular_velocity*dt;
        end
        
        % up and down
        if evalin('base','exist(''KEYS_SHIFT'',''var'') && KEYS_SHIFT')
            userCommands.velocityZ = userCommands.velocityZ + USER_PARAMS.acceleration*dt;
        elseif evalin('base','exist(''KEYS_CONTROL'',''var'') && KEYS_CONTROL')
            userCommands.velocityZ = userCommands.velocityZ - USER_PARAMS.acceleration*dt;
        else
            userCommands.velocityZ = 0;
        end

        userCommands.velocityZ = clip(userCommands.velocityZ, -5, 5);
        userCommands.pitch = clip(userCommands.pitch, -25*pi/180, 25*pi/180);
        userCommands.roll = clip(userCommands.roll, -25*pi/180, 25*pi/180);
        userCommands.psi = clip(wrapToPi(userCommands.psi),-0.995*pi,0.995*pi);
        
        % user commands are in body fixed frame
        phi = stateVector(1); theta = stateVector(2); psi = stateVector(3);        

        g = 9.81;       
        m = USER_PARAMS.mass_kg;           % kg 

        Kp = 0.75*m; 
        v_u  = stateVector(12);
        VelocityErrorLimit = 2.5;
        error_v_u = clip(userCommands.velocityZ - v_u, -VelocityErrorLimit, VelocityErrorLimit);
       
        userCommands.force = Kp*error_v_u + m*g; 

        force = userCommands.force;
        angleCommands = [userCommands.pitch; userCommands.roll; userCommands.psi];        
end

function setupKeyHandler(figHandle)
    % keyboard functionality determined via Chat GPT
    set(figHandle, 'KeyPressFcn',  @keyDown);
    set(figHandle, 'KeyReleaseFcn',@keyUp);
end

function keyDown(~, event)
    assignin('base', ['KEYS_' upper(event.Key)], true);
end

function keyUp(~, event)
    assignin('base', ['KEYS_' upper(event.Key)], false);
end

function startSimCB(fig)
    ui = getappdata(fig,'SIMDATA');
    if isempty(ui); ui = struct('running',false,'requestReset',false,'hFig',fig); end
    ui.running = true;
    ui.requestReset = false;
    setappdata(fig,'SIMDATA',ui);
end

function stopSimCB(fig)
    ui = getappdata(fig,'SIMDATA');
    if isempty(ui); ui = struct('running',false,'requestReset',false,'hFig',fig); end
    ui.running = false;
    setappdata(fig,'SIMDATA',ui);
end

function restartSimCB(fig)
    ui = getappdata(fig,'SIMDATA');
    if isempty(ui)
        ui = struct('running',false,'requestReset',false,'hFig',fig); 
    end
    ui.requestReset = true;     % main loop will re-init sim
    setappdata(fig,'SIMDATA',ui);
end

