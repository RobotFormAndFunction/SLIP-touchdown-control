function [out, tout]=interpmethod(robot, IC, sim_param,num,func_param,func_type,uncert_type)
% Runs the simulation with all the different interpolation methods
%data should be the optimal solution.  outputs how many succesful steps
%were taken
% Inputs:
%           - data: Data set found using optimization
%           - robot: robot class
%           - IC: are the system's initial conditions
%           - sim_param: system parameters
%           - Num: number of steps
%           - func_param: this is a structure with everything the lookup table needs
%           - func_type: [f1 f2] how you want things to work enter the
%               following for f1
%                 0) only perform the optmization
%                 1) only perform the interpolation
%                 2) perform both optimization and interpolation
%               for f2
%                 3) Add noise to optimization
%               - 0) performs lookup every step
%               - 1) performs lookup every x steps where x is defined
%                    by func_param.step.  This will initiate the count part
%                    in func_param
%               - 2) Does correction based on previous step data where how
%                    have back is determiend by func_param.step


if func_type(2)==0
    f=@td_table_lookup_var_ter;
elseif func_type(2)==1
    f=@td_skip;
elseif func_type(2)==2
    f=@td_old;
end

%run the simulation using optimization to find td angle
if func_type(1)==0 || func_type(1)==2
    func_param.type=2;
    data=Stance_sim_Ode(robot,IC,sim_param,num,0,f,func_param);

end

if func_type(1)==3
    
    data=[];
    func_param.type=3;
    
    for i=1:length(uncert_type)
        func_param.uncert=uncert_type(i);

        %finds td with interpolation

        data_lookup=Stance_sim_Ode(robot,IC,sim_param,num,0,f,func_param);
        data_lookup.tdcalc;

        %adds all the data together
        data=[data, data_lookup];

    end
end

if func_type(1)==1 || func_type(1)==2

    if func_type(1)==1
        data=[];
    end

    func_param.type=1;
    %methods={'linear' 'nearest', 'baseline'};
    methods={'linear', 'baseline'};

    for i=1:length(methods)
        func_param.typ=methods{i};

        if strcmp(methods{i},'baseline')
            %this is with no control
            sim_param.flightmodel=@flight_conservative_ode;
            data_lookup=Stance_sim_Ode(robot,IC,sim_param,num);
            h=floor(sum((data_lookup.switch==1))/2);
            data_lookup.tdcalc=ones(1,h)*IC(3);
            data_lookup.t   = [];
            data_lookup.b   = [];
            data_lookup.time= [];

        else

            %finds td with interpolation
            func_param.type=1;
            data_lookup=Stance_sim_Ode_p_v_t(robot,IC,sim_param,num,0,f,func_param);
            data_lookup.tdcalc; 
            
        end

        %adds all the data together
        data=[data, data_lookup];

        %plotting for paper writing
        %plot_stance(data,data_lookup,31)

    end
end

%compiles all the data together
[out,tout]=outputpad(num,data);

end

function [out,tout]=outputpad(num,x)
%outputs all matrix of td values for easy viewing x can be 1xM matrix of
%data structures from Stanc_Sim.  First row is how many steps were taken
%based on length of td calc.
l=length(x);

out=NaN(ceil(num)+2,l);
tout=NaN(ceil(num),l);

for i=1:l

    m=length(x(i).tdcalc);
    t=length(x(i).time);

    out(2:m+1,i)=x(i).tdcalc';

    tout(1:t,i)=x(i).time';

    out(1,i)=sum(x(i).switch)/2;

end

end

function plot_stance(data,data_lookup,num_steps)
figure(752)

subplot(2,1,1)
hold on
index=find(data.switch==1);
in=index(numsteps);
data.switch=data.switch(1:in);
plotvec=[1; find(data.switch==1)];
plotdata=[data.rec_states(:,1),data.rec_states(:,4)];

type=1;
for i=1:length(plotvec)-1
    if type ==1
        plot(plotdata(plotvec(i):plotvec(i+1),1),plotdata(plotvec(i):plotvec(i+1),2),'r', 'LineWidth', 3.5)
        plot(plotdata(plotvec(i):plotvec(i+1),1),ones(plotvec(i+1)-plotvec(i)+1,1)*data.yoffset(i),'k', 'LineWidth', 3.5)
        type=2;

    elseif type ==2
        plot(plotdata(plotvec(i):plotvec(i+1),1),plotdata(plotvec(i):plotvec(i+1),2),'b--', 'LineWidth', 3.5)
        plot(plotdata(plotvec(i):plotvec(i+1),1),ones(plotvec(i+1)-plotvec(i)+1,1)*data.yoffset(i),'k', 'LineWidth', 3.5)
        type=1;
    end
end
font_s=18;
xl = xlabel('X pos [m]');
yl = ylabel('Y pos [m]');
set(xl, 'FontSize', font_s)
set(yl, 'FontSize', font_s)
grid minor
ax = gca; % Get the current axes handle
ax.YColor='k';
set(ax, 'FontSize', font_s); % Adjust 14 to your desired font size

axis([-0.02 1.85 -0.01 0.04])
hold off


subplot(2,1,2)
hold on

index=find(data_lookup.switch==1);
in=index(numsteps);
data_lookup.switch=data_lookup.switch(1:in);
plotvec=[1; find(data_lookup.switch==1)];
plotdata=[data_lookup.rec_states(:,1),data_lookup.rec_states(:,4)];

type=1;
for i=1:length(plotvec)-1
    if type ==1
        plot(plotdata(plotvec(i):plotvec(i+1),1),plotdata(plotvec(i):plotvec(i+1),2),'r', 'LineWidth', 3.5)
        plot(plotdata(plotvec(i):plotvec(i+1),1),ones(plotvec(i+1)-plotvec(i)+1,1)*data_lookup.yoffset(i),'k', 'LineWidth', 3.5)
        type=2;

    elseif type ==2
        plot(plotdata(plotvec(i):plotvec(i+1),1),plotdata(plotvec(i):plotvec(i+1),2),'b--', 'LineWidth', 3.5)
        plot(plotdata(plotvec(i):plotvec(i+1),1),ones(plotvec(i+1)-plotvec(i)+1,1)*data_lookup.yoffset(i),'k', 'LineWidth', 3.5)
        type=1;
    end
end
%title('Point Mass Location')
font_s=18;
xl = xlabel('X pos [m]');
yl = ylabel('Y pos [m]');
set(xl, 'FontSize', font_s)
set(yl, 'FontSize', font_s)
grid minor
ax = gca; % Get the current axes handle
ax.YColor='k';
set(ax, 'FontSize', font_s); % Adjust 14 to your desired font size
lgd=legend('Stance', 'Terrain Height','Flight','Location','southoutside','Orientation','horizontal');
lgd.NumColumns = 3;
legend('boxoff')
fontsize(lgd,16,'points')
axis([-0.02 1.85 -0.01 0.04])
hold off


end
