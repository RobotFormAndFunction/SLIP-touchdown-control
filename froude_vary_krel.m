%For a purely conservative model tests various froude numbers
clear, close all, format compact

%% Parameter selection
mass=100;
l0=0.075;

krel_ob=linspace(10,25,6);


%Flags
save_flag=0;
make_data_set_flag=0;
data_centered=0; %set to 0 to run for entire data search area of td and beta.
%                 set to 1 to run a denser data set around the initial
%                 points

%System parameters
robot=Robot_Param;
robot.mass=mass/1000;
robot.l0=l0;
[sim_param, ~, ~] = make_sim_param('simulation parameters.txt');
sim_param.stancemodel=0;
sim_param.flightmodel=0;

fignum=1;

%% Make the data sets
if make_data_set_flag==1

    %create td and beta vectors
    load('m_1200g_krel_10_data_set.mat','magsol2','td_vec','beta_vec','IC');
    td_old=td_vec;
    beta_old=beta_vec;

    if data_centered == 1
        load('m_1200g_krel_10_IC.mat','IC');

        %Need to pass is beta angle (angle of velocity magnitude) and touchdown
        %angle
        t_range_mod = 0.175;
        b_range_mod = 0.175;

        trun = ceil((t_range_mod*2)/deg2rad(0.09));
        brun = ceil((b_range_mod*2)/deg2rad(0.09));

        %Find initial beta
        [~,dx,~,dy]=p2r(IC);
        b_init=atan2(dy,dx);


        %Create td and beta vectors
        td_vec    = linspace(IC(3)-t_range_mod,IC(3)+t_range_mod,trun);
        beta_vec  = linspace(b_init+b_range_mod,b_init-b_range_mod,brun);

    elseif data_centered==0
        nrun=350;
        td_vec    = linspace(0.4,pi/2-0.001,nrun);
        beta_vec  = linspace(-0.05,-pi/2+0.001,nrun);
    end

    for i=1:length(krel_ob)


        % Create the simulation parameters
        robot.k0=(krel_ob(i)*sim_param.g*(mass)/1000)/l0;

        % Perform optimization

        [tot_eng, tot_eng2, magsol, magsol2, step, step2]=data_set_centered_mex(robot,sim_param,td_vec,beta_vec);

        % Save data
        if save_flag==1
            if data_centered == 0
                data_set = strcat('m_',num2str(mass),'g_krel_',num2str(krel_ob(i)),'_krel_vary.mat');
            elseif data_centered==1
                data_set = strcat('m_',num2str(mass),'g_krel_',num2str(krel_ob(i)),'_krel_vary_cent.mat');
            end
            save(data_set);
        end



    end
end

%% Read in data otherwise

if make_data_set_flag==0

    [froude10,~,~] = load_froude_data(mass,10,data_centered);
    [froude13,~,~] = load_froude_data(mass,13,data_centered);
    [froude16,~,~] = load_froude_data(mass,16,data_centered);
    [froude19,~,~] = load_froude_data(mass,19,data_centered);
    [froude22,~,~] = load_froude_data(mass,22,data_centered);
    [froude25,beta_vec,td_vec] = load_froude_data(mass,25,data_centered);

end

%% Find Slopes only do if make_data_set_flag==0

if make_data_set_flag==0

    Ax_fr    = zeros(length(td_vec),length(beta_vec));
    bx_fr    = zeros(length(td_vec),length(beta_vec));

    Ax_cot    = zeros(length(td_vec),length(beta_vec));
    bx_cot    = zeros(length(td_vec),length(beta_vec));

    residuals = zeros(length(td_vec),length(beta_vec));

    for i=1:length(froude10(1,:))

        for j=1:length(froude10(:,1))

            y=[froude10(i,j), froude13(i,j), froude16(i,j), froude19(i,j), froude22(i,j), froude25(i,j)];
            [P,S]=polyfit(krel_ob,y,1);

            Ax_fr(i,j) = P(1);
            bx_fr(i,j) = P(2);
            residuals(i,j)=S.normr;

        end

    end

end

cplot(td_vec,-beta_vec,Ax_fr',linspace(0.001,2.5,100),1004,'Slope')
cplot(td_vec,-beta_vec,residuals',linspace(0.001,0.6,100),1005,'res')

%% Plot specific slices

i=120;
j=100;

i=40;
j=20;

i=77;
j=1;

y=[froude10(i,j), froude13(i,j), froude16(i,j), froude19(i,j), froude22(i,j), froude25(i,j)];
[P,S]=polyfit(krel_ob,y,1);
y2=polyval(P,krel_ob);
S.normr

figure(102)
close(102)
figure(102)
hold on
title('Froude')
plot(krel_ob,y,'r','LineWidth',3.5)
plot(krel_ob,y2,'k--','LineWidth',3.5)
xlabel('Krel'), ylabel('Fr')
hold off

%% Functions

function [froude, beta_vec,td_vec]=load_froude_data(mass,krel,centered)
%inputs mass is mass value
%       krel is the relative stiffness
%       centered should be 0 for non centered data creation or 1 for yes it
%       is centered

if centered == 0
    data_set = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_krel_vary.mat');
elseif centered==1
    data_set = strcat('m_',num2str(mass),'g_krel_',num2str(krel),'_krel_vary_cent.mat');
end

load(data_set,'magsol2', 'beta_vec', 'td_vec','robot');
froude=magsol2.^2/(9.81*robot.l0);

end

function cplot(x,y,z,c,fignum,title_var)
[x,y]=meshgrid(x,y);

figure(fignum)
hold on

contourf(x,y,z,c)
colorbar;
xlabel('Ground Angle \theta'), ylabel('Velocity Angle [$\beta$]','interpreter','latex')
title(title_var)
hold off

end

