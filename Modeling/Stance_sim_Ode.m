function [data]=Stance_sim_Ode(robot_param, IC,sim_param, iter, motor_parameter, func, func_param)
% STANCE_SIM will run a complete leg stance cycle starting from stance
% Inputs:
%   Robot_Param class
%
%   IC, a vector of initial conditions in the following form
%       - [zeta0, dzeta0, phi0, dphi0] with the following units
%       - [[m], [m/s], [rad], [rad/sec]]
%
%   iter, optional parameter, defaults to 1 if not provided, number of
%       stance/flight loops to run
%
%   simulation_param class
%
% motor_parameter a class only needed if a motor is being used in the simulation
%
%
% func, optional, third party function you want to call will take in:
%       n,data,simoutput,robot, IC, sim_param,motor_parameter,state
%       (flight, stance), func param
%       called during initilization and then after each simulation run
%
% func_param anything that needs to be passed into func, can be any data
%       type func needs to be set up to interpret
%
% Returns data structure with the following feilds:
%       - data.rec_states         vector of rec states [x,x',x",y,y',y"]
%       - data.polar_states       vector of polar states [phi, phi', phi", z, z', z"]
%                                   during flight all states are returned
%                                   as 0
%       - data.switch    vector of 0,1 1 denotes a switch between stance,
%                        flight
%       - data.yoffset   Shows nominal ground position
%       - data.PeGravity vector of potential energy in system from gravity
%       - data.PeSpring  vector of potential energy in system from the spring
%       - data.KE        vector of system kinetic energy
%       - data.Etot      total energy sum of PESpring, PEGravity, KE
%       - data.flightIC  vector of flight initial conditions
%       - data.stanceIC  vector of stance initial conditions
%   if motor ==1 the following will be returned
%       - data.voltage   motor voltage
%       - data.load_torq torque applied to motor
%       - data.omega     motor speed
%       - data.curr      motor current, may be 0 is LA not provided
%       - data.T_m       motor output torque

%Check robot_param structure and return values
[~, l0, ~, ~, xoffset, yoffset]=robot_param.unpack;

% Initial Conditions
if length(IC)~=4 error('Need 4 initial conditions\n'); end
sIC=IC; %rewrite to stance IC varaible to pervent loosing initial values

%Simulation Parameters
[~, ~, ~, stancemodel, flightmodel, ~, terrainvar, tervec]=sim_param.unpack;

%Set number of iterations
if ~exist('iter','var') iter=1;  end

% Motor Set Up
if exist('motor_parameter','var') & ~isa(motor_parameter, 'double')
    motor=1;
    [Ra, La, kt, kb, J, C, ~, ~, R]=motor_parameter.unpack;
    mvar=[Ra, La, kt, kb, J, C, R];
else
    motor_parameter=motor_param;
    motor=0;
    mvar=0;
    R=1;
end

if ~exist('func_param','var')
    func_param.inc=0;
end

%set up empty variables

data.rec_states   =  [];
data.polar_states =  [];
data.switch       =  [];
data.yoffset      =  [];
data.PeGravity    =  [];
data.PeSpring     =  [];
data.KE           =  [];
data.Etot         =  [];
data.flightIC     =  [];
data.stanceIC     =  IC;
data.phidot       =  [];
data.load_torq    =  [];
data.tvec         =  [];
if motor == 1
    if isequal(stancemodel,@stance_motor_ode) || isequal(stancemodel,@stance_motor_ode_poly) ...
            || isequal(stancemodel,@stance_motor_conV)
        sIC=[IC, motor_parameter.curr0];
    end
    data.voltage   = [];
    data.omega     = [];
    data.curr      = [];
    data.T_m       = []; %torque provided by the motor
end

%Call Third Party Function if it exists
if exist('func')
    if isa(func, 'function_handle')
        [data, func_param]=func(1,data,[],robot_param, IC, sim_param,motor_parameter, 0, func_param);
    end
end

if imag(IC(2))~=0
    dsvjkln=3;
end

modelType=1;

for i=1:2*iter
    
    switch modelType
        %% Stance Model
        case 1  %Run stance model
            
            optionsode=odeset('Abstol', 1e-8, 'Reltol', 1e-8, 'Events', @stance_switch);
            [tt,yout]=ode45(stancemodel,[0,2],sIC,optionsode,sIC,robot_param,sim_param, ...
                motor_parameter,func_param,mvar,motor_parameter.stanceV,motor_parameter.interpType);

            %Polar States
            yout_ac=stancemodel(tt,yout',sIC,robot_param,sim_param, motor_parameter,...
                func_param,mvar,motor_parameter.stanceV,motor_parameter.interpType);
            polar_states=[yout(:,3:4), yout_ac(4,:)', yout(:,1:2), yout_ac(2,:)'];
            
            %convert to rectangular states
            [x,dx,y,dy,ddx,ddy]=p2r(polar_states,1);
            rec_states=[xoffset+x,dx,ddx,yoffset+y,dy,ddy];
            if i > 2
                if (min(rec_states(:,4))<yoffset  ...
                        || abs(data.rec_states(end,4)-rec_states(1,4))>0.001)  && sim_param.opt==0
                    break
                end
            end
            
            %Create the exit flag
            exitflag=zeros(length(yout(:,1)),1); exitflag(end)=1;
            
            %Calculate Energy Terms
            KE=1/2*robot_param.mass*(rec_states(:,2).^2+rec_states(:,5).^2);
            PeGravity=robot_param.mass*sim_param.g*rec_states(:,4);
            PeSpring=1/2*robot_param.k0*(robot_param.l0-yout(:,1)).^2;
            Etot=KE+PeGravity+PeSpring;
            
            %set up IC may need to modify this if you write a new ODE func
            %as need to change the IC vec
            
            if isequal(flightmodel,@flight_conservative_ode)
                p=rec_states(end,:);
                fIC=[p(1),p(2),p(4),p(5), -IC(3)];
                
            elseif ~isequal(flightmodel,@flight_conservative_ode)
                p=rec_states(end,:);
                fIC=[p(1),p(2),p(4),p(5),2*pi-polar_states(end,1) ...
                    R*polar_states(end,2), motor_parameter.curr0];
            else
                p=rec_states(end,:);
                fIC=[p(1),p(2),p(4),p(5), -polar_states(1,end)];
            end
            if imag(fIC(1))~=0
                dsvjkln=3;
            end

            data.yoffset  = [data.yoffset; yoffset];
            
            %Handles variable terrain
            if terrainvar==1
                if tervec~=100
                    if i<2*iter-1
                        yoffset=tervec(i+1);
                    end
                else
                    yoffset=normrnd(0,l0/sim_param.yVar);
                    
                end
            end
            func_param.yoffset=yoffset;
            modelType=2;
            
        case 2
            %% Flight Model        
            optionsode=odeset('Abstol', 1e-8, 'Reltol', 1e-8, 'Events', @flight_switch);
            [tt,yout]=ode45(flightmodel,[0,2],fIC,optionsode,fIC,robot_param,sim_param,...
                motor_parameter,func_param,mvar,motor_parameter.flightV,motor_parameter.interpType);
            
            %Get acc values
            yout_ac=flightmodel(tt,yout',fIC,robot_param,sim_param, motor_parameter,...
                func_param,mvar,motor_parameter.flightV,motor_parameter.interpType);
            rec_states=[yout(:,1:2), yout_ac(2,:)', yout(:,3:4), yout_ac(4,:)'];
            if i>3
                if (rec_states(end,4)<data.yoffset(end)  ...
                        || abs(data.rec_states(end,4)-rec_states(1,4))>0.001)  && sim_param.opt==0
                    break
                end
            end
            
            KE=1/2*robot_param.mass*(rec_states(:,2).^2+rec_states(:,5).^2);
            PeGravity=robot_param.mass*sim_param.g*rec_states(:,4);
            PeSpring=zeros(length(rec_states(:,1)),1);
            Etot=KE+PeGravity;
            
            %Create the exit flag
            exitflag=zeros(length(yout(:,1)),1); exitflag(end)=1;
            
            %Set IC for stancemodel
            
            phi=-yout(end,5);
            
            p=yout(end,1:4);
            x    = -robot_param.l0*cos(phi); 
            xvel = p(2); 
            y    = robot_param.l0*sin(phi); 
            yvel = p(4);
            dzeta0 = (x*xvel+y*yvel)/sqrt(x^2+y^2);
            dphi0  = -(x*yvel-y*xvel)/(x^2+y^2);
            sIC     = [robot_param.l0,dzeta0,phi,dphi0];
            xoffset       = rec_states(end,1)-x;
            data.yoffset  = [data.yoffset; yoffset];
            
             if isequal(stancemodel,@stance_motor_ode) || isequal(stancemodel,@stance_motor_ode_poly) ...
                     || isequal(stancemodel,@stance_motor_conV)
                sIC=[sIC, motor_parameter.curr0];
             end
            
            %Write Polar states as 0;
            polar_states      =  zeros(length(rec_states(:,1)),6);
            polar_states(:,1) = yout(:,5);
            
            if length(yout(1,:))==7
                polar_states(:,2)=yout(:,6)/R;
            end
            
            modelType=1;
            
    end
   
    
    %% Format variable to output data structure
    
    %Save the rec and polar states
    data.rec_states =  [data.rec_states; rec_states];
    data.polar_states =  [data.polar_states; polar_states];
    
    %Save the time vector
    if i==1
        data.tvec=tt(:,1);
    else
        data.tvec=[data.tvec; data.tvec(end)+tt(:,1)];
    end
    
    %Switch
    data.switch =  [data.switch; exitflag];
    
    %Energy Data
    data.PeGravity  =  [data.PeGravity; PeGravity];
    data.PeSpring   =  [data.PeSpring; PeSpring];
    data.KE         =  [data.KE; KE];
    data.Etot       =  [data.Etot; Etot];
    
    % Other misc a few things only appear in one or the other
    if modelType==2
        data.flightIC = [data.flightIC; fIC];
        if motor == 0
            %Note does not account for a gearbox otherwise would be
            %m*L^2/gearR
            load_torq=(robot_param.mass*polar_states(:,1)).*(sim_param.g*cos(polar_states(:,4))...
                -2*polar_states(:,5).*polar_states(:,2));
        else
            load_torq=(robot_param.mass*polar_states(:,4).^2).*(2*polar_states(:,5).*polar_states(:,4)+...
                sim_param.g*cos(polar_states(:,1))./polar_states(:,4));
        end
       
        data.load_torq  = [data.load_torq; load_torq(:,1)];
    else
        data.stanceIC = [data.stanceIC; sIC(1:4)];
        temp=zeros(length(rec_states(:,5)),1);
        data.phidot = [data.phidot; temp];
        data.load_torq  = [data.load_torq; temp]; %toruq on motor currently assumed to be zero in flight
    end
    if motor == 1
        if modelType==2

            if isequal(stancemodel,@stance_motor_ode_poly)

                V=motor_parameter.stanceV;

                volts=zeros(length(tt),1);
                for z=1:length(tt)
                    for j=1:length(V)
                        volts(z)=volts(z)+V(j)*tt(z)^(j-1);
                    end
                end

                data.voltage=[data.voltage; volts];

            elseif isequal(stancemodel,@stance_motor_conV)

                data.voltage = [data.voltage; motor_parameter.stanceV(1)*ones(length(tt),1)];

            else

                data.voltage = [data.voltage; ...
                interpn(motor_parameter.stanceV(:,1),motor_parameter.stanceV(:,2),tt,motor_parameter.interpType)];

            end
            

        else
            if isequal(flightmodel,@flight_motor_ode_con_v) 

                data.voltage = [data.voltage; ones(length(tt),1)*motor_parameter.flightV(1)];
            
            elseif isequal(flightmodel,@flight_motor_ode_full_speed_for_t)
                
                volts=(tt<motor_parameter.flightV(1))*motor_parameter.flightV(2);
                data.voltage=[data.voltage; volts];

            else

                data.voltage = [data.voltage; ...
                    interpn(motor_parameter.flightV(:,1),motor_parameter.flightV(:,2),tt,motor_parameter.interpType)];

            end
        end
        
        data.omega = [data.omega; polar_states(:,2)*motor_parameter.gearR]; %motor shaft speed
        data.curr  = [data.curr; yout(:,end)];  %current should be last column of output
        data.T_m   = [data.T_m; yout(:,end)*motor_parameter.gearR*motor_parameter.kt]; %torque provided by the motor = curr *kt*gearR
        
        %Set new motor IC's
        motor_parameter.curr0  = data.curr(end);
        
    end

    %% Call Third Party Function if it exists
    if exist('func')
        if isa(func, 'function_handle')
            [data, func_param]=func(2,data,[],robot_param, IC, sim_param,motor_parameter,modelType, func_param);
        end
    end
    

end
end

%% Functions

function [value, isterminal, direction]=flight_switch(t,x,IC,robot,simp, motor_parameter,func_param,mvar,V,type)
if x(4)<0

    value = x(3)-robot.l0*sin(-x(5)) - func_param.yoffset;
    isterminal = 1;   % Stop the integration
    direction  = -1;   %z needs to be increasing
else
    value=1;
    isterminal=0;
    direction=0;
    
end
end

function [value, isterminal, direction]=stance_switch(t,x,IC,robot,simp, motor_parameter,func_param,mvar,V,type)
value      = (x(1) >= robot.l0);
isterminal = 1;   % Stop the integration
direction  = 1;   %z needs to be increasing
end