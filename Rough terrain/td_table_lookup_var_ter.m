function [data, func_param]=td_table_lookup_var_ter(n,data,simoutput,robot, IC, sim_param,motor_parameter,state, func_param)
%Table lookup to solve for TD angle on rough terrain requires the following
%fields in func_param:
%       type   - selects the type of calc you want to do:
%                   1) does a table lookup to solve for td at each step
%                   2) runs an optimization for step by step td angle
%                   3) Applies a perturbation to the optimal result
%       xvel   - x velocity of COM as this should be constant
%       beta   - vector of beta angles that corresponds to tabeng
%       td     - vector of touchdown angles that corresponds to tabeng
%       V0     - matrix of velocities for initial conditions
%       up_var - [optional] allows you to select how often energy is
%       updated specifically x velocity maybe you can only sample onece
%       ever 5 steps?


%% Write td variable to phi0
if n==1
    %sets up initial conditions for simulation
    data.tdcalc=[];
    func_param.cycle=0;
    data.t=[];
    data.b=[];
    data.time=[];
    return
end


%% Interpolation
if func_param.type==1 && state==2

    % if func_param.correctxvel==1
    %     xvel_TO=data.rec_states(end,2);
    % else
    %     xvel_TO=func_param.xvel;
    % end

    %grab data needed to find touchdown angle
    magsol   = func_param.V0; % pull out the initial conditions matrix
    td_vec   = func_param.td;
    beta_vec = func_param.beta;

    %Initial varaible declarations

    t_start=tic;

    %determine takeoff conditions
    xvel_TO  = data.rec_states(end,2);
    yvel_TO  = data.rec_states(end,5);
    y_TO     = data.rec_states(end,4);
    TO_angle = data.polar_states(end,1); %might not be needed
    vmag_TO  = sqrt(yvel_TO^2+xvel_TO^2); %velocity magnituede at takeoff

    if xvel_TO<=0
        feasible=0;
    else
        feasible=1;
    end

    if feasible==1

        % Determine apex conditions of flight phase
        t_apex=yvel_TO/sim_param.g;            %time to get to max height
        y_apex=yvel_TO^2/(2*sim_param.g)+y_TO; %max height
        x_apex=t_apex*xvel_TO;                 %x pos of apex

        %Find ground distribution
        del_y=y_apex-robot.l0; %max ground height
        ground=linspace(del_y,-3*robot.l0*1/func_param.yVar,11); %ground vector

        %ground=data.yoffset(end);

        mag_interp=zeros(1,length(td_vec));
        td_lookup=zeros(1,length(ground));

        %%

        for jj=2:length(ground)
        %for jj=1:length(ground)

            td_heights=ground(jj)+robot.l0*sin(td_vec);

            deltaPE_vec=2*9.81*(y_TO-td_heights);

            v_vec=sqrt(vmag_TO^2+deltaPE_vec);
            yvel_vec=-sqrt(v_vec.^2-xvel_TO.^2);
            beta_desire=atan2(yvel_vec,xvel_TO);

            for k=1:length(beta_desire)
                mag_interp(k)=interp1(beta_vec,magsol(k,:),beta_desire(k));
            end

            mag_err=(mag_interp-v_vec).^2./v_vec;

            % make sure we don't accidentally catch one of the data points
            % on the transition from infeasible to feasible
            a=findpeaks(mag_err);
            if isempty(a) || a(1)>length(td_vec)-10
                in=1;
            else
                in=find(mag_err==a(1));
            end

            %find index that minimizes the error
            inl=find(mag_err(in:end)==min(mag_err(in:end)));

            if isempty(inl)
                inl=1;
            end

            td_lookup(jj)=td_vec(inl+in-1);

        end

        td_lookup(1)=td_lookup(2);
        td=mean(td_lookup);

        %find theta vs t
        ypos_vec = y_apex-ground(2:end)-robot.l0*sin(td_lookup(2:end));

        t_vec = sqrt(ypos_vec/(0.5*sim_param.g))+t_apex; %
        [P]=polyfit(t_vec,td_lookup(2:end),2);


        % plotting for the purpose of paper writing
        plot_f=0;
        if plot_f==1
            tp=polyval(P,t_vec);
            figure()
            hold on
            plot(t_vec,td_lookup(2:end), "r*",'MarkerSize',10)
            plot(t_vec,tp,'k','LineWidth',3.5)

            font_s=18;
            xl = xlabel('Time [s]');
            yl = ylabel('\theta_{TD} [rad]');

            grid minor

            set(xl, 'FontSize', font_s)
            set(yl, 'FontSize', font_s)
            ax = gca; % Get the current axes handle
            set(ax, 'FontSize', font_s); % Adjust 14 to your desired font size

            lgd = legend("Solved points","Curve Fit");
            legend('boxoff')
            fontsize(lgd,14,'points')
            hold off
        end

        func_param.apex_t=t_apex;
        func_param.td_start=td_lookup(1);
        func_param.p_v_t=P;

    else
        td=data.stanceIC(end,3);
    end



    t_end=toc(t_start);


end

%% Optimization
if (func_param.type==2 || func_param.type==3) && state==2

    num_iter=40;

    %func_param.yoff=data.yoffset(end);

    sim_param_new=sim_param;

    %make copy of robot param and assign start height to last end height
    robot_new=robot;
    robot_new.yoffset=data.yoffset(end);

    sim_param.terrainvar=1;
    sim_param_new.tervec=[data.yoffset(end), evalin('caller','yoffset')*ones(1,2*num_iter-1)];


    t_start=tic;
    [td, ~] = fminbnd(@td_opt,0.4,1.5,[],robot_new,sim_param_new,data.stanceIC(end,:),3);

    if func_param.type==3
        td=td+normrnd(0,func_param.uncert);
    end

    t_end=toc(t_start);
    % debugging
    %func_param.yoff=evalin('caller','yoffset');
    %func_param.td=td;
    %data_1=Stance_sim_Ode(robot_new,ICnew,sim_param_new,2,0,@con_td_write,func_param);
    %fprintf('Steps taken at optimal: %f \n', sum(data_1.switch)/2)


    % tvec=linspace(0.4,ub,200);
    % fv=zeros(1,200);
    % step=zeros(1,200);
    % for z=1:length(tvec)
    %     fv(z)=td_opt(tvec(z),robot_new,sim_param_new, ICnew,5,func_param);
    %     func_param.td=tvec(z);
    %     data_td=Stance_sim_Ode(robot_new,ICnew,sim_param_new,num_iter,0,@con_td_write,func_param);
    %     %graphing(data_td,sim_param_new,3,1,[data_td.rec_states(:,1),data_td.rec_states(:,4)], ["r", "b--", "k"])
    %     step(z)=sum(data_td.switch)/2;
    %     if step(z)==num_iter
    %         break
    %     end
    % end
    % index=find(step==max(step));
    % td=tvec(index(1));
    %
    % % figure
    % % hold on
    % % axis([0, pi/2, -40, 10])
    % % ylabel('Cost Function Value')
    % % xlabel('Touchdown Angle [Rad]')
    % % plot(tvec,fv)
    % % yyaxis right
    % % plot(tvec,step)
    % % ylabel('Steps to Failure')
    % % hold off
    % fprintf('Steps taken at optimal: %f %f  TD angles %f, %f \n', [sum(data_1.switch)/2, step(index(1)) td1 td])
    % data_1.switch(1)=1;
    % index=find(data_1.switch==1);
    % ty=1;
    % figure
    % hold on
    % for h=1:length(index)-1
    %     if ty ==1
    %         ty=2;
    %         plot(data_1.rec_states(index(h):index(h+1),1), ...
    %         data_1.rec_states(index(h):index(h+1),4), 'r','LineWidth',1.5)
    %     else
    %         ty=1;
    %         plot(data_1.rec_states(index(h):index(h+1),1), ...
    %         data_1.rec_states(index(h):index(h+1),4), 'b--','LineWidth',1.5)
    %     end
    %
    % end
    % hold off

end



if state==2
    % Write Output to caller function
    data.tdcalc=[data.tdcalc, td];

    data.time=[data.time, t_end];

    tempIC=evalin('caller','fIC');
    tempIC(5)=-td;

    assignin('caller','fIC',tempIC)
end

end
