%For a purely conservative model what is the energy of the system that
%would result in a stable gait for a givcen touch down angle and velocity
%angle
clc, clear, close all, format compact

%Perturbing a stable model and then finding a new touchdown angle to
%correct for the perturbation

clc, clear, close all, format compact

%%Parameter selection
mass = 50;
krel = 22;
l0   = 0.02;

fignum=1000;

%% Flags
IC_flag               = 1; %set to 1 to run an optimization for IC otherwise 0 reads in old IC
td_variability_flag   = 0; %flag for finding variability in td angle
rough_opt_flag        = 1; %flag for performing opt pertabation
rough_opt_w_noise_flag= 0; %flag for performing opt pertabation with noise on touchdown angle
rough_LU_flag         = 0; %flag for performing LU pertabation
generate_terrain_flag = 0; %flag for generating new terrain
save_flag             = 1; %saves data
interp_data_flag      = 1; %perform interpolation or read in existing data set


%% Create the simulation parameters and read in data
if interp_data_flag==0
    data_set_filename        = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_data_set.mat');
end

%Male simulaton parameters
[sim_param, ~, ~] = make_sim_param('simulation parameters.txt');

robot=Robot_Param;
robot.mass = mass/1000;
robot.l0   = l0;
robot.k0   = (krel*robot.mass*9.81)/l0;
robot.bl   = 0;

sim_param.stancemodel=@stance_model;
sim_param.flightmodel=@flight_conservative_ode;

ter_var_filename  = strcat('m_',num2str(mass),'g_l0_',num2str(robot.l0*1000),'mm_terrain.mat');

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

%% Rough Terrain setup
% Setup vectors

if interp_data_flag==0
    [beta_vec,td_vec,~,~,magsol,step,beta_step,td_step]=load_td_data(data_set_filename,2,2); %this is needed just to populate the functon call
elseif interp_data_flag==1
    [beta_vec,td_vec,magsol,beta_step,td_step]=extrapolate_system_energy(deg2rad(0.09*64),deg2rad(0.09*64),krel,robot);
end


% Reset values for rough terrain
[x,dx,y,dy] = p2r(IC);
beta        = atan2(dy,dx);
mag         = sqrt(dx^2+dy^2);

%setup simulation parameters
rough_val.V0=magsol; %initial velocity matrix
rough_val.td=td_vec;
rough_val.beta=beta_vec;
rough_val.xvel=dx;
rough_val.type=2;
rough_val.typ='spline';
rough_val.step=1; %how many steps to take between values

sim_param.terrainvar=1;
sim_param.yVar=15;
c=10;

%% Find TD variability on rough terrain
%this basically varies the ground roughness and then sees how much td needs
%to change to reamin xstable on that ground. Gives a sense of how much the
%robot needs to correct the td angle and gives an idea of stability
if td_variability_flag == 1
    fprintf('Running TD_var for M=%2d, Krel=%2d\n' ,[mass,krel])

    yVar=1./[0.2 0.18 .16 .14 .12 .1 .08 0.07 .06 0.055 0.05 0.045 .04 0.035 .03 0.025 .02];

    %find new IC for optimization
    sim_param.terrainvar=0;
    data=Stance_sim_Ode(robot,IC,sim_param,td_var_num,0,@td_table_lookup_var_ter,rough_val);
    sim_param.terrainvar=1;
    graphing(data,sim_param,3,fignum,[data.rec_states(:,1),data.rec_states(:,4)], ["r", "b--", "k"])
    fignum=fignum+1;
    IC=data.stanceIC(end,:);

    td_std=zeros(length(yVar),1);

    for k=1:length(yVar)
        sim_param.yVar=yVar(k);
        data=Stance_sim_Ode(robot,IC,sim_param,td_var_num,0,@td_table_lookup_var_ter,rough_val);
        m=mean(data.tdcalc(1:sum(data.switch)/2-1));
        td_std(k)=std(data.tdcalc(1:sum(data.switch)/2-1)-m);
        fprintf('yVar: %f Steps: %f\n', [yVar(k), sum(data.switch)/2])

    end
    fprintf('----------------------------------\n')

    %remove NaN values
    for w=1:length(td_std)
        if isnan(td_std(w))
            td_std(w)=min(td_std);
        end
    end

    figure(fignum)
    hold on
    plot(1./yVar,rad2deg(td_std),'*')
    xlabel('Terrain Height As a Fraction of Leg Length')
    ylabel('Std of Td angle needed for stability [Deg]')
    axis([0.02,0.12, 0, 4])
    hold off
    fignum=fignum+1;

    if save_flag==1
        td_var_filename = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_td_var.mat');
        save(td_var_filename,'td_std','yVar')
    end


end

%% Looking at different interpolation methods
%This part uses 8 different methods of table lookups and returns
% how many steps the robot was able to take for each of those conditions.
% It also records which method was able to take the most successful steps
% and tabulates that.  Finally it shows how many steps the random
% pertabation was behind the hightest value.  You can adjust the number of
% iterations with the Z variable and the number of simulated steps with
% num, and what yVar values you want to change with yVar.  Would like to
% functionalize but would be trickeyx
num=50; % number of simulated  steps
z=100;    % how many iterations to run of num steps change back to 100

yVar=1./[0.25 0.2 0.18 .16 .14 .12 .1 .08 0.07 .06 0.055 0.05 0.045 .04 0.035 .03 0.025 .02];

if generate_terrain_flag==1
    terrain=create_terrains(robot,100,100,yVar);
else
    load(ter_var_filename,'terrain','yVar')
end

%performs optimization for rough terrain
if rough_opt_flag == 1
    fprintf('Running Opt Rough terrain for M=%2d, Krel=%2d\n' ,[mass,krel])
    dim=0;
    typ=0;
    cycle_type=0;

    rough_val.std=zeros(length(yVar),1);

    %Run everything
    [av_step,av_t]=rough_interp(cycle_type, typ, num, z, yVar, terrain,...
        robot, IC, sim_param,rough_val,0);

    if save_flag==1
        rough_opt_filename = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_Rough_opt.mat');
        save(rough_opt_filename,'av_step','av_t')
    end

end

if rough_opt_w_noise_flag == 1
    fprintf('Running Opt with noise Rough terrain for M=%2d, Krel=%2d\n' ,[mass,krel])
    typ=3;
    cycle_type=0;

    sigma=deg2rad(0.09*[1 2 4 8 16 32 64]);

    rough_val.std=zeros(length(yVar),1);

    %Run everything
    [av_step,av_t]=rough_interp(cycle_type, typ, num, z, yVar, terrain,...
        robot, IC, sim_param,rough_val,sigma);

    if save_flag==1
        rough_opt_filename = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_noise_Rough_opt.mat');
        save(rough_opt_filename,'av_step','av_t')
    end

end

%performs LU for rough terrain
if rough_LU_flag == 1
    fprintf('Running LU Rough terrain for M=%2d, Krel=%2d\n' ,[mass,krel])

    dim=1;
    typ=1;
    cycle_type=[0 0];

    %Run everything
    [av_step,av_t]=rough_interp(cycle_type, typ, num, z, yVar, terrain,...
        robot, IC, sim_param,rough_val,0);

    if save_flag==1
        td_var_filename = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_Rough_LU.mat');
        save(td_var_filename,'av_step','av_t')
    end


end


%% Functions

function [av_step,av_t]=rough_interp(cycle_type, typ, num, z, yVar, terrain,...
    robot, IC, sim_param,rough_val,uncert_val)

c=0;

%typ=0 for optimal and 1 for lookup 3 for opt with uncertainity
if typ==0
    dim=1;
elseif typ==1
    dim=2;
elseif typ==2
    dim=3;
elseif typ==3
    dim=length(uncert_val);
end

steptot = zeros(z,dim,length(yVar),length(cycle_type));
stepout = zeros(z*(num+1),dim,length(yVar),length(cycle_type));
tout    = NaN(z*num,dim,length(yVar),length(cycle_type));
t_here=tic;

for i=1:length(yVar)
    rough_val.yVar=yVar(i);

    for k=1:z

        %This creates a terrain.  So that all methods are tested on the
        %same terrain
        sim_param.tervec=terrain(k,:,i);

        for j=1:length(cycle_type)
            if j>1 %changed this
                rough_val.correctxvel=1;
            else
                rough_val.correctxvel=0;
            end

            %graphing(data,sim_param,3,1023,[data.rec_states(:,1),data.rec_states(:,4)], ["r", "b--", "k"])
    
            [b1, t1] = interpmethod(robot, IC, sim_param,num,...
                rough_val,[typ, cycle_type(j)],uncert_val); %everystep
            steptot(k,:,i,j) = b1(1,:);
            stepout((k-1)*(1+num)+1 : k*(1+num),:,i,j) = b1(2:end,:);
            tout((k-1)*num+1 : k*num,:,i,j) = t1;
            c=c+1;

        end
    end
    fprintf("Done with %d of %d\n", [i, length(yVar)])

end

% Split up data
av_step = zeros(2,dim,length(yVar),length(cycle_type));
av_t    = zeros(2,dim,length(yVar),length(cycle_type));

for j=1:length(cycle_type)
    for i=1:length(steptot(1,1,:,j))

        av_step(1,:,i,j) = mean(steptot(:,:,i,j));
        av_t(1,:,i,j)    = mean(tout(:,:,i,j),"omitnan");

        av_step(2,:,i,j) = std(steptot(:,:,i,j));
        av_t(2,:,i,j)    = std(tout(:,:,i,j),"omitmissing");

    end
end

toc(t_here)
end


function terrain=create_terrains(robot,ternum, num,yVar)
%creates a ternum number of terrains that are num long
% yVar corresponds to how many layers deep the matrix is in the third dimension
terrain=[zeros(ternum,num*2,length(yVar))];

for i=1:length(yVar)
    for j=1:length(terrain(:,1,1))
        for zizzles=3:2:2*num-1
            yoff=normrnd(0,robot.l0/yVar(i));
            terrain(j,zizzles:zizzles+1,i)=yoff;
        end
    end
end

end






