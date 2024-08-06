function [suc_total, av_step, av_time]=td_pertabation_V2(robot,sim_param,IC,iterations,rand_par,step_suc, type, beta_vec,td_vec,magsol)
% Startiing with and IC applies a perturbation to both velocity magnitude
% and direction. System takes one step with perturbation then recalculates
% the TD angle for step 2 that tries to give a stable gait.
% Solves for a new Td angle that returns a stable gait does this using 3
% different methods. Returns the number of iterations that we able to
% complete given number of steps.
%
% Solves using 4 methods:
%      1) Finds the td angle using optimization by trying to keep x' and y'
%         constant between the start and end of stance + the number of
%         steps taken
%
%      2) Assumes velocity magnitude is known at takeoff
%
%      3) BASELINE does not correct the TD angle and just holds constant
%
%Inputs:
%       - robot: Robot parameter class
%       - sim_param:  Simulation parameter class
%       - IC:         Initial conditions
%       - iterations: How many cycles you want it to run
%       - rand_par:   [vector] Sigma of the random pertabation applied to the system
%                     is mag/rand_per(i), or beta/rand_par(i) so a larger
%                     number correlates to greater varaition.  Can be a
%                     vector of values.
%       - step_suc:   Minimum number of steps to be considered a success
%       - type:       Methods used, enter 0 for [1,2], 1 for [1], or 2
%                     for [3].  Allows only optimization, all three or
%                     only lookup to be performed. 3 for [4]
%       - beta_vec:   velocity beta angle vector
%       - td_vec:     Td vector
%       - magsol:     Velocity matrix
%
% Outputs:
%      - Suc_total: matrix where first column is 1/rand_par(i).  Columns
%                   2,3,4 are are the the number of iterations that made it
%                   over step_suc.  Last column is number of unrecoverable
%                   trials
%      - av_step:   Columns 1,2,3 are the average number of
%                   steps completed.  4,5,6 is the std_deviation of steps
%                   completed


c=10;

print_flag=0;


if type==0
    if print_flag==1
        fprintf('    Opt, Interp,\n')
    end
    suc=[0 0 0]; %how many were greater than 10 steps, last counts number of unrecoverable
    av_step   = zeros(length(rand_par),4);
    av_time   = zeros(length(rand_par),2);
    suc_total = zeros(length(rand_par),4);
    step_vec=zeros(iterations,2);
    t_vec=zeros(iterations,2);

else
    if print_flag==1
        if type==1
            fprintf('    Opt\n')
        elseif type==2
            fprintf('    Baseline\n')
        else
            fprintf('    Baseline\n')
        end
    end
    suc=[0 0]; %how many were greater than 10 steps second number is how many unrecoverable steps there were
    av_step   = zeros(length(rand_par),2);
    av_time   = zeros(length(rand_par),1);
    suc_total = zeros(length(rand_par),3);
    step_vec=zeros(iterations,1);
    t_vec=zeros(iterations,1);

end

outs=1;

for j=1:length(rand_par)
    rnd_per=rand_par(j); % rnd_per is the sigma level applied where sigma is nominal/rnd_per
    for i=1:iterations

        %solve for velocity magnitude and direction
        [x,dx,y,dy]=p2r([robot.l0,IC(2),IC(3),IC(4)]);
        beta=atan2(dy,dx);
        mag=sqrt(dx^2+dy^2);

        %apply a random pertebation
        mag_new=normrnd(mag,abs(mag/rnd_per));
        beta_new=normrnd(beta,abs(beta/rnd_per));

        %create new IC's from pertabations
        xvel=mag_new*cos(beta_new);
        yvel=mag_new*sin(beta_new);

        ICnew=IC;
        ICnew(2)=(x*xvel+y*yvel)/sqrt(x^2+y^2);
        ICnew(4)=-(x*yvel-y*xvel)/(x^2+y^2);

        % Run one stance cycle to find exit criteria
        data=Stance_sim_Ode(robot,ICnew,sim_param,0.5);

        td_heights=robot.l0*sin(td_vec); %find heights

        %Grab takeoff conditions
        TO_angle = data.polar_states(end,1);
        xvel_TO=data.rec_states(end,2);
        yvel_TO=data.rec_states(end,5);
        vmag_TO=sqrt(xvel_TO^2+yvel_TO^2);

        %Check if system can clear step
        if (yvel_TO^2-2*9.81*(robot.l0-robot.l0*sin(TO_angle)))<0 || xvel_TO<0
            recoverable=0;
            td_lookup=IC(3);
            td=IC(3);
            data_lookup=data;
            opt_t=0;
            lk_t=0;
            suc(end)=suc(end)+1;
        else
            recoverable=1;
        end

        if (type==0 || type==2) && recoverable==1

            % start calcs specific to calculation
            tic

            deltaPE_vec=2*9.81*(robot.l0*sin(TO_angle)-td_heights);

            v_vec=sqrt(vmag_TO^2+deltaPE_vec);
            yvel_vec=-sqrt(v_vec.^2-xvel_TO.^2); %need to take the negative because the square doesn't consider the sign


            beta_desire=atan2(yvel_vec,xvel_TO);
            mag_interp=zeros(1,length(beta_desire));

            for k=1:length(beta_desire)
                mag_interp(k)=interp1(beta_vec,magsol(k,:),beta_desire(k));
            end


            %Plot outputs for paper writing

            plot_flag=0;

            if plot_flag==1
                figure(86)
                hold on
                plot(td_vec,beta_desire,"LineWidth",3.5)
                font_s = 18;
                xl = xlabel('Touchdown Angle \theta_{TD} [rad]');
                yl = ylabel('Velocity Angle \beta [rad]');
                set(xl, 'FontSize', font_s)
                set(yl, 'FontSize', font_s)
                ax = gca; % Get the current axes handle
                set(ax, 'FontSize', font_s); % Adjust 14 to your desired font size
                hold off

                figure(87)
                hold on
                plot(td_vec,v_vec,"LineWidth",3.5)
                font_s = 18;
                xl = xlabel('Touchdown Angle \theta_{TD} [rad]');
                yl = ylabel('Velocity magnitude [m/s]');
                set(xl, 'FontSize', font_s)
                set(yl, 'FontSize', font_s)
                ax = gca; % Get the current axes handle
                set(ax, 'FontSize', font_s); % Adjust 14 to your desired font size
                plot(td_vec,mag_interp,'ro',"MarkerSize",3.5,'MarkerFaceColor','r')
                lgd = legend('Theoretical','Data set');
                legend('boxoff')
                fontsize(lgd,14,'points')
                hold off
            end


            %Find velocity error and touchdown angle that corresponds to
            %minumum energy
            mag_err=(mag_interp-v_vec).^2./v_vec;
            td_lookup=(td_vec(mag_err==min(mag_err)));

            lk_t=toc;
        end

        %% Run Optimization
        if (type==0 || type==1) && recoverable == 1
            %Finish setting up some flight conditions
            robot.bl=0;
            sim_param.terrainvar=0;
            sim_param.stancemodel=@stance_model;
            sim_param.flightmodel=@flight_conservative_ode;

            %Need to pass a beta angle (angle of velocity magnitude) and
            %touchdown angle
            tic
            [td, ~]=fminbnd(@td_opt,0.4,1.5,[],robot,sim_param,ICnew,3);
            opt_t=toc;

            % Change plot_cost value from 0 to 1 to perform an exhaustive search
            % and plot outputs

            plot_cost=0;

            if plot_cost==1

                tvec=linspace(0.1,pi/2,250);
                fv=zeros(1,250);
                step=zeros(1,250);
                for z=1:length(tvec)
                    fv(z)=td_opt(tvec(z),robot,sim_param, ICnew,3);
                    func_param.td=tvec(z);
                    data=Stance_sim_Ode(robot,ICnew,sim_param,30,0,@con_td_write,func_param);
                    step(z)=sum(data.switch)/2;
                end

                figure(c*100)
                hold on
                axis([0, pi/2, 0, 10])
                ylabel('Cost Function Value')
                xlabel('Touchdown Angle [Rad]')
                plot(tvec,fv)
                yyaxis right
                plot(tvec,step)
                ylabel('Steps to Failure')
                hold off

            end
        end

        % Run_simulation

        if (type==0 || type==1 || type==3) && recoverable==1
            if type==3
                td=IC(3);
            end
            func_param.td=td;
            data= Stance_sim_Ode(robot,ICnew,sim_param,30,0,@con_td_write,func_param);

        end

        if (type==0 || type==2) && recoverable==1

            func_param.td=td_lookup;
            data_lookup=Stance_sim_Ode(robot,ICnew,sim_param,30,0,@con_td_write,func_param);

        end

        c=c+1;

        if type==0

            p=[td td_lookup sum(data.switch)/2 sum(data_lookup.switch)/2 ];
            if p(3)>step_suc suc(1)=suc(1)+1; end
            if p(4)>step_suc suc(2)=suc(2)+1; end

            step_vec(i,:)=[p(3), p(4)];
            t_vec(i,:)=[opt_t,lk_t];

        else
            if type==1 || type==3
                p=[td sum(data.switch)/2];

                if type==3
                    opt_t=0;
                end

                t_vec(i)=opt_t;

            else
                p=[td_lookup sum(data_lookup.switch)/2];
                t_vec(i)=lk_t;
            end


            if p(2)>step_suc suc(1)=suc(1)+1; end

            step_vec(i)=p(2);


        end


    end

    %Format outputs and print data to terminal

    if print_flag==1

        fprintf('----------------------------------------------------\n')

        if type==0

            fprintf('Sigma value: %1.5f \n', 1/rnd_per)
            fprintf('Number of steps>5 %1.5f %1.5f nonrecoverable steps %1.5f\n',suc)
            suc_total(outs,:)=[1/rnd_per, suc];
            suc= [0 0];

        else

            fprintf('Sigma value: %1.5f \n', 1/rnd_per)
            fprintf('Number of steps>5  %1.5f nonrecoverable steps %1.5f\n',suc)
            suc_total(outs,:)=[1/rnd_per, suc];
            suc= 0;

        end

        fprintf('----------------------------------------------------\n')

    else

        if type==0
            suc_total(outs,:)=[1/rnd_per, suc];
            suc= [0 0 0];
        else
            suc_total(outs,:)=[1/rnd_per, suc];
            suc= [0 0];
        end

    end

    av_step(outs,:)=[mean(step_vec), std(step_vec)];
    av_time(outs,:)=mean(t_vec);
    outs=outs+1;

end
end