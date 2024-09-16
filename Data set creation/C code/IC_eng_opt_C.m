function [aof]=IC_eng_opt_C(mag,robot_param,sim_param, beta, td, opt_type)
% IC_eng_opt
%Given a velocity angle at touchdown, and touchdown angle solves for a
%system energy value that results in stable gait. 
%
%INPUTS:
%   - robot_param : robot parameter class
%   - sim_param   : simulation parameter class
%   - beta        : velocity angle as defined as positive in the clockwise
%                   rotation with zero being the axis parrallel to the
%                   horizon that passes through the point mass
%   - td          : touchdown angle
% opt_type 1 optimizes the energy change between start and end
%          2 minimizes the change IC's (dphi0, zeta0)
%          3 min error between xdot and y dot at start and end of stance


%Use the magnitude, beta, and theta values and convert that into zdot and
%phidot values
[x,~,y,~]=p2r_C([robot_param.l0,0,td,0],0);

xvel=mag*cos(beta);
yvel=mag*sin(beta);

%Convert from x,y vel to polar
zdot=(x*xvel+y*yvel)/sqrt(x^2+y^2);
phidot=-(x*yvel-y*xvel)/(x^2+y^2);

stancemodel = @stance_model_C;
flightmodel = @flight_conservative_ode_C;

func_param.yoffset=0;

IC=[robot_param.l0,zdot,td,phidot];

switch opt_type
    
    case 1
        data=Stance_sim_Ode_C(robot_param,IC,sim_param,1,stancemodel,flightmodel);
        if min(data.rec_states(end,4)) < 0
            aof=abs(data.Etot(1)-data.Etot(end))*10000;
        else
            aof=abs(data.Etot(1)-data.Etot(end));
        end
        
    case 2
        data=Stance_sim_Ode_C(robot_param,IC,sim_param,1,stancemodel,flightmodel);
        if min(data.rec_states(end,4)) < 0
            aof=sum(abs(IC-data.stanceIC))*1000;
        else
            aof=abs(IC(2)-data.stanceIC(2))+abs(IC(4)-data.stanceIC(4));
        end
        
    case 3

        sim_param.opt=0;
        data=Stance_sim_Ode_C(robot_param,IC,sim_param,1,stancemodel,flightmodel);

        %for a 0.5, 1.5 etc number of steps need to convert the second
        %part from - to +

        if isempty(data.polar_states)
            aof=1000;
        else
            index=find(data.switch==1);
            aof=1000*(data.rec_states(1,2)-data.rec_states(index(end),2))^2+1000*(data.rec_states(1,5)-data.rec_states(index(end),5))^2;
            if min(data.rec_states(end,4))< 0 || data.switch(end)~=1
                aof=aof*1000;
            end
        end

    case 4

        sim_param.opt=0;
        data=Stance_sim_Ode_C(robot_param,IC,sim_param,1,stancemodel,flightmodel);

        %for a 0.5, 1.5 etc number of steps need to convert the second
        %part from - to +

        aof=10*(data.rec_states(1,2)-data.rec_states(end,2))^2+10*(data.rec_states(1,5)-data.rec_states(end,5))^2;
    case 5

        sim_param.opt=1;
        data=Stance_sim_Ode_C(robot_param,IC,sim_param,1,stancemodel,flightmodel);

        %for a 0.5, 1.5 etc number of steps need to convert the second
        %part from - to +

        aof=10*(data.rec_states(1,2)-data.rec_states(end,2))^2+10*(data.rec_states(1,5)-data.rec_states(end,5))^2;

        
end

end