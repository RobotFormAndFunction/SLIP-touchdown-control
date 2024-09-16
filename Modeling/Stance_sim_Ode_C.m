function [data]=Stance_sim_Ode_C(robot_param, IC,sim_param, iter, stancemodel, flightmodel)
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

sIC=IC; %rewrite to stance IC varaible to pervent loosing initial values

%Simulation Parameters
[~, ~, ~, ~, ~, ~, terrainvar, tervec]=sim_param.unpack;

motor=0;

%set up empty variables

dataLen=20000*2*iter;

data.rec_states   =  zeros(dataLen,6);
data.polar_states =  zeros(dataLen,6);
data.switch       =  zeros(dataLen,1);
data.yoffset      =  zeros(2*iter,1);
data.PeGravity    =  zeros(dataLen,1);
data.PeSpring     =  zeros(dataLen,1);
data.KE           =  zeros(dataLen,1);
data.Etot         =  zeros(dataLen,1);
data.flightIC     =  zeros(2*iter,5);
data.stanceIC     =  zeros(2*iter,4);
data.stanceIC(1,:)=  IC;
data.phidot       =  zeros(dataLen,1);
data.load_torq    =  zeros(dataLen,1);
data.tvec         =  zeros(dataLen,1);
fIC               =  [0, 0, 0, 0, 0];

datLen=0;

rec_states=zeros(20000,6);
polar_states=zeros(20000,6);

KE=zeros(20000,1);
PeGravity=zeros(20000,1);
PeSpring=zeros(20000,1);
Etot=zeros(20000,1);

tt=zeros(20000,1);
yout=zeros(20000,4);
yout1=zeros(20000,5);

dataStorageLoc=1;
dataStorageLoc_old=1;
fcount=0;
scount=0;

coder.varsize('rec_states')

modelType=1;

for i=1:2*iter
    
    switch modelType
        %% Stance Model
        case 1  %Run stance model
            
            optionsode=odeset('Abstol', 1e-8, 'Reltol', 1e-8, 'Events', @stance_switch);
            [tt,yout]=ode45(stancemodel,[0,2],sIC,optionsode,sIC,robot_param,sim_param);

            %Polar States
            yout_ac=stancemodel(tt,yout',sIC,robot_param,sim_param);
            polar_states=[yout(:,3:4), yout_ac(4,:)', yout(:,1:2), yout_ac(2,:)'];
            
            datLen=length(yout(:,1));

            %convert to rectangular states
            [x,dx,y,dy,ddx,ddy]=p2r_C(polar_states,1);
            rec_states=[xoffset+x,dx,ddx,yoffset+y,dy,ddy];
            if i > 2
                if (min(rec_states(:,4))<yoffset  ...
                        || abs(data.rec_states(dataStorageLoc_old,4)-rec_states(1,4))>0.001)  && sim_param.opt==0
                    break
                end
            end
            
            %Calculate Energy Terms
            KE=1/2*robot_param.mass*(rec_states(:,2).^2+rec_states(:,5).^2);
            PeGravity=robot_param.mass*sim_param.g*rec_states(:,4);
            PeSpring=1/2*robot_param.k0*(robot_param.l0-yout(:,1)).^2;
            Etot=KE+PeGravity+PeSpring;
            
            %set up IC may need to modify this if you write a new ODE func
            %as need to change the IC vec
            
            if isequal(flightmodel,@flight_conservative_ode_C)
                p=rec_states(end,:);
                fIC=[p(1),p(2),p(4),p(5), -IC(3)];
                
            else
                p=rec_states(end,:);
                fIC=[p(1),p(2),p(4),p(5), -polar_states(1,end)];
            end

            data.yoffset(i)  = yoffset;
            
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
            modelType=2;
            fcount=fcount+1;
            
            
        case 2
            %% Flight Model        
            optionsode=odeset('Abstol', 1e-8, 'Reltol', 1e-8, 'Events', @flight_switch);
            [tt,yout1]=ode45(flightmodel,[0,2],fIC,optionsode,fIC,robot_param,sim_param,yoffset);
            
            %Get acc values
            yout_ac=flightmodel(tt,yout1',fIC,robot_param,sim_param,yoffset);
            rec_states=[yout1(:,1:2), yout_ac(2,:)', yout1(:,3:4), yout_ac(4,:)'];
            
            datLen=length(yout1(:,1));

            if i>3
                if (rec_states(end,4)<data.yoffset(i)  ...
                        || abs(data.rec_states(dataStorageLoc_old,4)-rec_states(1,4))>0.001)  && sim_param.opt==0
                    break
                end
            end
            
            KE=1/2*robot_param.mass*(rec_states(:,2).^2+rec_states(:,5).^2);
            PeGravity=robot_param.mass*sim_param.g*rec_states(:,4);
            PeSpring=zeros(length(rec_states(:,1)),1);
            Etot=KE+PeGravity;
            
            %Set IC for stancemodel
            
            phi=-yout1(end,5);
            
            p=yout1(end,1:4);
            x    = -robot_param.l0*cos(phi); 
            xvel = p(2); 
            y    = robot_param.l0*sin(phi); 
            yvel = p(4);
            dzeta0 = (x*xvel+y*yvel)/sqrt(x^2+y^2);
            dphi0  = -(x*yvel-y*xvel)/(x^2+y^2);
            sIC     = [robot_param.l0,dzeta0,phi,dphi0];
            xoffset       = rec_states(end,1)-x;
            data.yoffset(i)  = yoffset;
            
            
            %Write Polar states as 0;
            polar_states      =  zeros(length(rec_states(:,1)),6);
            polar_states(:,1) = yout1(:,5);
            
            modelType=1;
            scount=scount+1;
            
    end
   
    
    %% Format variable to output data structure
    
    %Save the rec and polar states
    datpEnd=datLen+dataStorageLoc-1;
    data.rec_states(dataStorageLoc:datpEnd,:)   = rec_states;
    data.polar_states(dataStorageLoc:datpEnd,:) = polar_states;
    
    %Save the time vector
    if i==1
        data.tvec(dataStorageLoc:datpEnd)=tt(:,1);
    else
        data.tvec(dataStorageLoc:datpEnd,:)=data.tvec(dataStorageLoc)+tt(:,1);
    end
    
    %Switch
    data.switch(datpEnd) = 1;
    
    %Energy Data
    data.PeGravity(dataStorageLoc:datpEnd) =  PeGravity;
    data.PeSpring(dataStorageLoc:datpEnd)  =  PeSpring;
    data.KE(dataStorageLoc:datpEnd)        =  KE;
    data.Etot(dataStorageLoc:datpEnd)      =  Etot;
    
    % Other misc a few things only appear in one or the other
    if modelType==2
        data.flightIC(fcount,:) = fIC;
        if motor == 0
            %Note does not account for a gearbox otherwise would be
            %m*L^2/gearR
            load_torq=(robot_param.mass*polar_states(:,1)).*(sim_param.g*cos(polar_states(:,4))...
                -2*polar_states(:,5).*polar_states(:,2));
        else
            load_torq=(robot_param.mass*polar_states(:,4).^2).*(2*polar_states(:,5).*polar_states(:,4)+...
                sim_param.g*cos(polar_states(:,1))./polar_states(:,4));
        end
       
        data.load_torq(dataStorageLoc:datpEnd)  = load_torq;
    else
        data.stanceIC(scount,:) = sIC(1:4);
    end

    dataStorageLoc_old=dataStorageLoc;
    dataStorageLoc=datpEnd+1;
end
end

%% Functions

function [value, isterminal, direction]=flight_switch(t,x,IC,robot,simp,yoffset)
if x(4)<0

    value = x(3)-robot.l0*sin(-x(5)) - yoffset;
    isterminal = 1;   % Stop the integration
    direction  = -1;   %z needs to be increasing
else
    value=1;
    isterminal=0;
    direction=0;
    
end
end

function [value, isterminal, direction]=stance_switch(t,x,IC,robot,simp)
value      = (x(1) >= robot.l0);
isterminal = 1;   % Stop the integration
direction  = 1;   %z needs to be increasing
end