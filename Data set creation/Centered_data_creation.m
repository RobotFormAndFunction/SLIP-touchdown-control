%For a purely conservative model creates a really fine mesh of data points around the
%system's initial condition
clc, clear, close all, format compact

%% Parameter selection
mass_vec=[120];
krel_vec=[15];

%Flags
save_flag=1;

%initial varaible selection
options=optimoptions('fmincon','Display','off','Algorithm','sqp','OptimalityTolerance',1e-8,'ConstraintTolerance',1e-8,...
    'StepTolerance',1e-8);

fignum=1;
count=1;
en=4;
system_num=1;

%%

for ji=1:length(mass_vec)
    mass=mass_vec(ji);
    for jj=1:length(krel_vec)
        krel=krel_vec(jj);

        %% Create the simulation parameters and read in data

        data_set_filename        = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_data_set.mat');
        IC_filename              = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_IC.mat');
        pertabation_LU_filename  = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_LU_results.mat');
        ter_var_filename         = strcat('m_',num2str(mass),'g_terrain.mat');

        [sim_param, ~, ~] = make_sim_param('simulation parameters.txt');
        
        %Load in IC's, magnitude, beta,td vectors

        load(data_set_filename,'robot');

        load(IC_filename,'IC');

        % check initial condition
        data=Stance_sim_Ode(robot,IC,sim_param,100);
        fprintf('Number of Steps IC takes: %f \n', sum(data.switch)/2)
        graphing(data,sim_param,3,fignum,[data.rec_states(:,1),data.rec_states(:,4)], ["r", "b--", "k"])
        fignum=fignum+1;
        if sum(data.switch)/2<90
            error('IC not sufficient')
        end

        %% Perform optimization

        % Range calculations
        t_range_mod = 0.175;
        b_range_mod = 0.175;

        trun = ceil((t_range_mod*2)/deg2rad(0.09));
        brun = ceil((b_range_mod*2)/deg2rad(0.09));

        t1=tic();
        % trun=5;
        % brun=5;

        %Find initial beta
        [~,dx,~,dy]=p2r(IC);
        b_init=atan2(dy,dx);

        %Create td and beta vectors
        td_vec    = linspace(IC(3)-t_range_mod,IC(3)+t_range_mod,trun);
        beta_vec  = linspace(b_init+b_range_mod,b_init-b_range_mod,brun);

        % td_vec    = linspace(1.1,1.4,trun);
        % beta_vec  = linspace(-0.4,-0.7,brun);

        td_vec(td_vec >0.4 & td_vec<pi/2);
        beta_vec(beta_vec >-pi/2 & beta_vec>-0.05);

        fprintf("Total number of elements %d\n", length(td_vec)*length(beta_vec))
        
        %run optimization
        sim_param.flightmodel=0;
        sim_param.stancemodel=0;
        [tot_eng, tot_eng2, magsol, magsol2, step, step2]=data_set_centered_mex(robot,sim_param,td_vec,beta_vec);

        %% Save data
        if save_flag==1
            data_set = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_data_set_fine.mat');
            save(data_set);
        end

        toc(t1)

        %% Plotting
        c=linspace(0,max(max(tot_eng)),25);
        cplot(td_vec,-beta_vec,tot_eng2',c,fignum+1,'Total Energy')


        c=linspace(0,max(max(magsol2)),25);
        cplot(td_vec,-beta_vec,magsol2',c,fignum+2,'Velocity Magnitude')

        c=linspace(1.5,100,50);
        cplot(td_vec,-beta_vec,step2',c,fignum+4,'Steps Taken')

        system_num=system_num+1;

        

    end
end

%% Functions

function cplot(x,y,z,c,fignum,title_var)
[x,y]=meshgrid(x,y);

figure(fignum)
hold on

contourf(x,y,z,c)
colorbar;
xlabel('Ground Angle [$\theta$]','interpreter','latex'), ylabel('Velocity Angle [$\beta$]','interpreter','latex')
title(title_var)
hold off

end

function [data, IC]=run_sim(robot,sim_param,vel_mag,td,beta,nrun)
%Convert coordinates

[x,~,y,~]=p2r([robot.l0,0,td,0]);

xvel=vel_mag*cos(beta);
yvel=vel_mag*sin(beta);

%Convert from x,y vel to polar
zdot=(x*xvel+y*yvel)/sqrt(x^2+y^2);
phidot=-(x*yvel-y*xvel)/(x^2+y^2);

IC=[robot.l0,zdot,td,phidot];

data=Stance_sim_Ode(robot,IC,sim_param,nrun);

end
