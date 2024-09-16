function [tot_eng, tot_eng2, magsol, magsol2, step, step2]=data_set_centered(robot, sim_param,td_vec,beta_vec)

%% Parameter selection

%initial varaible selection
options=optimoptions('fmincon','Display','off','Algorithm','sqp','OptimalityTolerance',1e-8,'ConstraintTolerance',1e-8,...
    'StepTolerance',1e-8);
count=1;
en=4;

%% Perform optimization

% lower and upper bounds for energy magnitude
lb=0;
ub=10;

fprintf("Total number of elements %d\n", int32(length(td_vec)*length(beta_vec)));

magsol    = zeros(length(td_vec),length(beta_vec));
IC_vec    = zeros(length(td_vec)*length(beta_vec),4);
tot_eng   = zeros(length(td_vec),length(beta_vec));
%per_ken   = zeros(length(td_vec),length(beta_vec));
step      = zeros(length(td_vec),length(beta_vec));

for i=1:length(td_vec)
    td=td_vec(i);
    for j=1:length(beta_vec)
        beta=beta_vec(j);

        % run optimization
        [magout, fval]=fmincon(@(mag)IC_eng_opt_C(mag,robot,sim_param,beta,td,3),en,[],[],[],[],lb,ub,[],options);

        % rum simulation with final results
        [data, IC]=run_sim(robot,sim_param,magout,td,beta,30);

        if sum(data.switch)<3
            %try other IC's if completes less than 3 steps
            ig=[0.5,1,1.5,2.5];
            for jim=1:length(ig)
                [mag1, fval1]=fmincon(@(mag)IC_eng_opt_C(mag,robot,sim_param,beta,td,3),ig(jim),[],[],[],[],lb,ub,[],options);
                [data1, IC1]=run_sim(robot,sim_param,mag1,td,beta,30);

                if sum(data1.switch)>sum(data.switch)
                    data=data1;
                    IC=IC1;
                    magout=mag1;
                    fval=fval1;
                end

            end

        end


        %Save data
        step(i,j)=sum(data.switch)/2;
        magsol(i,j)=magout;
        IC_vec(count,:)=IC;

        %Get the total system energy
        [x,~,y,~]=p2r_C(IC,0);
        tot_eng(i,j) = 1/2*robot.mass*magout^2+sim_param.g*robot.mass*y;

        % eng=1;
        % 
        % %find percent of kinetic energy at center of stance
        % if eng==1
        %     data=Stance_sim_Ode_C(robot,IC,sim_param,0.5,stancemodel,flightmodel);
        %     if ~isempty(data.polar_states)
        %         pos=find(data.polar_states(:,4)==min(data.polar_states(:,4)));
        %         per_ken(i,j) = data.KE(pos)/(data.Etot(1)-data.PeGravity(pos))*100;
        %     end
        % end
        count=count+1;

        if mod(count,100)==0
            fprintf('Done with %6d', int32(count));
            fprintf(' of %d\n',int32(length(td_vec)*length(beta_vec)));
        end
    end
end

%Remove data from any system that did not complete 2 steps
tot_eng2=tot_eng; tot_eng2 (step<1.5) = 0;
magsol2=magsol;   magsol2  (step<1.5) = 0;
step2=step;       step2    (step<1.5) = 0;

%remove any point with crazy high velocity (not sure I ever use these)
magsol3=magsol;   magsol3 (magsol>6) = 0;
tot_eng3=tot_eng; tot_eng2(magsol>6) = 0;

end


function [data, IC]=run_sim(robot,sim_param,vel_mag,td,beta,nrun)
%Convert coordinates
stancemodel = @stance_model_C;
flightmodel = @flight_conservative_ode_C;
[x,~,y,~]=p2r_C([robot.l0,0,td,0],0);

xvel=vel_mag*cos(beta);
yvel=vel_mag*sin(beta);

%Convert from x,y vel to polar
zdot=(x*xvel+y*yvel)/sqrt(x^2+y^2);
phidot=-(x*yvel-y*xvel)/(x^2+y^2);

IC=[robot.l0,zdot,td,phidot];

data=Stance_sim_Ode_C(robot,IC,sim_param,nrun,stancemodel,flightmodel);

end
