% DESCRIPTION:
% This script tests the implementation of the ALVA LET model for inferring
% layer moduli of spheroid cell clusters based from a optical 
% interferometry MicroPiette Aspiration (MPA) test. 

% The model consideres two layeres measured, where the top layer is
% measured from expermients.

% Layer E-modulus and the rigid body displacement are treated as unknowns.
% The rigid body displacements help guiding the optimization
% algorithm and is therfore included 

clear all, clc

% -------------------------------------------------------------------------
% Select response analysis type
% -------------------------------------------------------------------------
alva.analysis = 'Full';
% 1) 'Full'     : Conventional full integration with one-step Richardson
%                 extrapolation (for improved convergence near surface)
%                 applied for evaluation of displacements and stresses
% 2) 'PolFit'   : Use polynomial fit technique to reduce integral length
%                 according to Andersen et al. (2018) for evaluation of
%                 surface displacements

% -------------------------------------------------------------------------
% Select interface
% -------------------------------------------------------------------------
alva.bond = 'Bonded';
% 1) 'Bonded'       : Full bonding between layers
% 2) 'Slip'         : Interface bonding factor
% 3) 'Frictionless' : No bonding between layers

% -------------------------------------------------------------------------
% Numerical parameters
% -------------------------------------------------------------------------
alva.N  = 300;  % Number of Bessel zero points in numerical integration
alva.n  = 30;   % Number of Gauss points points between zero points.

%--------------------------------------------------------------------------
% Define path and filename for experimental data
%--------------------------------------------------------------------------
profile_path = '/Users/massimilianoberardi/Desktop/Displ_profile/T24/14';
profile_folders = dir(profile_path);

Eres = zeros((length(profile_folders)-3)*2,3);
E_idx = 1;

for k = 4:length(profile_folders)
    
    alva.path     = strcat(profile_path,'/',profile_folders(k).name);
    alva.filename = 'P';
    j=15;
    for i = 1:1%2
        %--------------------------------------------------------------------------
        % Import experimental data
        %--------------------------------------------------------------------------
        try
        exp = load(append(alva.path,'/',alva.filename,'_',num2str(j),'.mat'));
        j = 10 + 5;
        alva.expno = i;
        % -------------------------------------------------------------------------
        % Load configuration
        % -------------------------------------------------------------------------
        alva.a  = [26                     % Load radii [mm] (probe nozzle inner egde)
                   28.5                   % Load radii [mm] (probe nozzle inner egde)
                   26]*1e-3;              % Load radii [mm] (probe nozzle circumference)
    
        A1 = pi*(alva.a(1))^2;            % Load area  1 & 3 [mm2]
        A2 = pi*(alva.a(2))^2-A1;         % Load area  2 [mm2]
    
        alva.q = zeros(length(alva.a),1); % Initialize load vector
    
        alva.q(1)  = -exp.load*1e-6;      % Suction pressure [MPa] (probe nozzle)
        alva.q(2)  = -alva.q(1)*(A1/A2);  % Positive pressure [MPa] (probe nozzle circumference)
        alva.q(3)  = -alva.q(2);          % Negative pressure [MPa] (remove positive part for probe nozzle area)
    
        alva.Xl = [0.0 0.0
            0.0 0.0
            0.0 0.0];                     % Load positions [mm]: [x1 y1; x2 y2;..xi yi];
    
        % -------------------------------------------------------------------------
        % Sensor locations along z-axis in [mm]
        % -------------------------------------------------------------------------
    
        alva.fopos{i} = exp.zpos*1e-3; % from micrometer to millimeter
        
        % -------------------------------------------------------------------------
        % Response measurements (vertical surface displacements, uz) in [mm]
        % -------------------------------------------------------------------------
    
        alva.dz_exp{i} = exp.zdisp*1e-6; % from nanometer to millimeter
    
        % -------------------------------------------------------------------------
        % Spheroid size
        % -------------------------------------------------------------------------
    
        alva.spr{i} = exp.diameter*1e-3*0.5;  % Spheroid size/radius from micrometer to millimeter
        
        % -------------------------------------------------------------------------
        % Rigid body displacement
        % -------------------------------------------------------------------------
    
        alva.rbd{i} = min(alva.dz_exp{i});
    
        % -------------------------------------------------------------------------
        % Initialize optimization problem
        % -------------------------------------------------------------------------
        % Stop criteria
        tolfun  = 0.1e-6;             % Object function stop criteria
        tolvar  = 0.1e-6;             % Step size value stop criteria
        tolfval = 2e3;              % Maximum function evaluations
        tolit   = 6e3;              % Maximum iterations
        
        options = optimset('MaxIter',tolit,...
        'MaxFunEvals',tolfval,'TolFun',tolfun,'TolX',tolvar); %'Display','iter',

    % -------------------------------------------------------------------------
    % Spheroid material properties (minimum two layers required)
    % -------------------------------------------------------------------------
    nlf     = 1;                         % number of finite layers
    tmeas   = 20e-3;                     % measured layer thicknesss ~ 25 micron
    tz      = [tmeas];%ones(1,nlf)*tmeas;         % layer thicknesses of finite layers 
    alva.zi = cumsum(tz);                % Depth of first n-1 layers from the surface [mm]
    nz      = length(alva.zi)+1;         % number of layers

    alva.E  = ones(1,nz)*1e-3;           % Layer Young's moduli [MPa]   
    alva.nu = ones(1,nz)*0.5;           % Layer Poisson's ratio [-]
    alva.kh = ones(1,nz-1)*1e9;          % Interface bonding/horizontal spring [MPa/mm]
    alva.nz = nz;
    alva.tz = tz;

    % Initial parameters
    E0     = alva.E;        
    rbd0   = alva.rbd{i};
    tz0    = tz;
            
    % Constraints
    E_LB    = E0.*0.1;          
    E_UB    = E0.*20; 
    rbd_LB  = rbd0*0.97;
    rbd_UB  = rbd0*1.03;
    tz_LB   = ones(1,nz-1)*10e-3;      % min. thickness of one cell ~ 6 micron
    tz_UB   = tz_LB*5;                 
    
    % Transform search parameters (for increased sped/efficiency)
    x0    = [log10(E0*1e6), rbd0*1e3, tz0*1e2]; 
    x0_LB = [log10(E_LB*1e6), rbd_LB*1e3, tz_LB*1e2];   % Lower bound values
    x0_UB = [log10(E_UB*1e6), rbd_UB*1e3, tz_UB*1e2];   % Upper bound values
    
    tic
    % Run optimization algorithm
    [xmin,fval,exitflag,output] = fminsearchbnd(@(x)inv_loop_mpa_cons(x,alva),x0,x0_LB,x0_UB,options); %,options
    toc
    
    fprintf('Iterations %d: ',output.iterations)
    fprintf('Relative error %f\n',fval)
    
        % -------------------------------------------------------------------------
        % Get response for optimimal predicted E-moduli
        % -------------------------------------------------------------------------
    alva.Xd = [alva.fopos{i}.*0 alva.fopos{i}.*0 alva.fopos{i}];
    alva.E  = 10.^(xmin(1:nz))./1e6;  % transform to real value in [MPa]
    alva.rbd_inf = xmin(nz+1)*1e-3;   % transform to real value in [mm] 
    tz       = xmin(nz+2:end)*1e-2;
    alva.zi  = cumsum(tz);

    alva     = init_LET(alva);
        
        UZ(:,i)  = alva.uz - alva.rbd_inf; 
        Z(:,i)   = alva.Xd(:,3);
        E1(i)    = alva.E(1)*1e6;
        E2(i)    = alva.E(2)*1e6;
        RBD(i)   = alva.rbd_inf*1e3;
        R2fit = 1-sum(((alva.dz_exp{1}+UZ(:,1))).^2)/sum(((alva.dz_exp{1}+mean(UZ(:,1)))).^2)

        Eres(E_idx,1:2) = alva.E;
        Eres(E_idx,3) = R2fit;
        E_idx = E_idx + 1
        catch
        Eres(E_idx,1:2) = alva.E;
        Eres(E_idx,3) = R2fit;
        E_idx = E_idx + 1
        end
    end
    %%
    figure(2), title('Displacement profile');
hold on, grid on
plot(alva.dz_exp{1}*1e3,alva.fopos{1}*1e3,'ko','LineWidth',1.25,'MarkerSize',5)
hold on
plot(-UZ(:,1)*1e3,Z(:,1)*1e3,'-.k','LineWidth',1.0,'MarkerSize',6)
   
end


%%
Eres_good = Eres(Eres(:,3)>0.75,1:2)