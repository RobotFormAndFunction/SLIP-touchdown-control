%Perturbing a stable model and then finding a new touchdown angle to
%correct for the perturbation

clc, clear, close all, format compact

%%Parameter selection
mass = 50;
krel = 22;
l0   = 0.02;

num=250; %number of trials to run
fignum=1000;

%% Flags
IC_flag          = 1; %set to 1 to run an optimization for IC otherwise 0 reads in old IC
pert_LU_flag     = 1; %flag for performing LU pertabation
opt_flag         = 0; %flag for performing opt pertabation
baseline_flag    = 0; %flag for performing baseline pertabation
save_flag        = 1; %flag to save data
interp_data_flag = 1; %if zero uses exsiting data set or interpretes data set if set to 1


%% Create the simulation parameters and read in data
if interp_data_flag==0
    data_set_filename        = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_data_set.mat');
end
pertabation_LU_filename  = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_LU_results.mat');

%Male simulaton parameters
[sim_param, ~, ~] = make_sim_param('simulation parameters.txt');

robot=Robot_Param;
robot.mass = mass/1000;
robot.l0   = l0;
robot.k0   = (krel*robot.mass*9.81)/l0;
robot.bl   = 0;

sim_param.stancemodel=@stance_model;
sim_param.flightmodel=@flight_conservative_ode;

%% Find an Initial Condition
options=optimoptions('fmincon','Display','off','TolFun',1e-8);

if IC_flag==1
    % [-1 9]
    IC=[robot.l0  	-3	1	0];
    [d,s]=fmincon(@SLIPopt,[-1.7 25],[],[],[],[],[-10 ,0],[0, 100],[],options, ...
        robot,sim_param,IC(3),3);
    IC=[robot.l0, d(1), IC(3), d(2)];
end

data=Stance_sim_Ode(robot,IC,sim_param,100);
index=find(data.switch==1);
fprintf('Number of Steps IC takes: %f \n', sum(data.switch)/2)
graphing(data,sim_param,3,fignum,[data.rec_states(:,1),data.rec_states(:,4)], ["r", "b--", "k"])
fignum=fignum+1;
if sum(data.switch)/2<25
    error('IC not sufficient')
end

%% Calculate system energy
[~,dx,~,dy]=p2r(IC);
beta=atan2(dy,dx);
fprintf('dx=%2.6f, dy=%2.6f, Beta %2.6f\n', [dx dy beta])

%% Perterb input and solve for new TD angle to retain stability

if interp_data_flag==0
    [beta_vec,td_vec,~,~,magsol,step,beta_step,td_step]=load_td_data(data_set_filename,2,2); %this is needed just to populate the functon call
elseif interp_data_flag==1
    [beta_vec,td_vec,magsol,beta_step,td_step]=extrapolate_system_energy(1,1,12,robot);
end


error_vec=1./[0.2 0.18 .16 .14 .12 .1 .08 0.07 .06 0.055 0.05 0.045 .04 0.035 .03 0.025 .02];

%this performs all of the optimization outputs
if opt_flag == 1
    t_start=tic;
    [step_suc_opt, av_step_opt,av_time_opt]=td_pertabation_V2(robot, sim_param,IC,num,error_vec, 6, 1, beta_vec,td_vec,magsol);
    toc(t_start)
    fprintf('M=%2d, Krel=%2d Opt done \n', [mass, krel])

    if save_flag==1
        pertabation_opt_filename = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_pert_opt.mat');
        save(pertabation_opt_filename);
    end
end

%this performs all of the baseline calculations
if baseline_flag == 1
    t_start=tic;
    [step_suc_base, av_step_base,av_time_base]=td_pertabation_V2(robot, sim_param,IC,num,error_vec, 6, 3, beta_vec,td_vec,magsol);
    toc(t_start)
    fprintf('M=%2d, Krel=%2d Baseline done \n', [mass, krel])

    if save_flag==1
        pertabation_baseline = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_baseline.mat');
        save(pertabation_baseline);
    end
end

% go through and do all of the interpolation
if pert_LU_flag == 1

    if interp_data_flag==0
        t_in=[1,2,3,5,6,10, 15, 25, 30, 50, 75];
    elseif interp_data_flag==1
        t_in=deg2rad(0.09)*[1 2 4 8 16 32 64];
        %t_in=deg2rad(0.09)*[2];
    end


    %t_in=2;
    L=length(error_vec);
    step_suc    = zeros(L*length(t_in),3);
    av_step     = zeros(L*length(t_in),2);
    av_tim      = zeros(L*length(t_in),1);
    td_step_vec = zeros(L*length(t_in),1);
    b_step_vec  = zeros(L*length(t_in),1);


    % Perform Lookup
    t_start=tic;
    for k=1:length(t_in)

        if interp_data_flag==0
            [beta_vec,td_vec,~,~,magsol,~,beta_step,td_step]=load_td_data(data_set_filename,t_in(k),t_in(k));
        elseif interp_data_flag==1
            [beta_vec,td_vec,magsol,beta_step,td_step]=extrapolate_system_energy(t_in(k),t_in(k),krel,robot);
        end

        [step_t, av_step_t,av_time_t]=td_pertabation_V2(robot, sim_param,IC,num,error_vec, 6, 2, beta_vec,td_vec,magsol);

        b_step_vec (L*(k-1)+1:k*L,:)  = beta_step;
        td_step_vec(L*(k-1)+1:k*L,:)  = td_step;
        step_suc   (L*(k-1)+1:k*L,:)  = step_t;
        av_step    (L*(k-1)+1:k*L,:)  = av_step_t;
        av_tim     (L*(k-1)+1:k*L,:)  = av_time_t;
        fprintf('Done with %2.0f of %2.0f \n', [k,length(t_in)])

    end
    toc(t_start)

    if save_flag==1
        save(pertabation_LU_filename);
    end
end








