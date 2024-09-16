function [beta_vec_new,td_vec_new,V0,beta_step,td_step]=extrapolate_system_energy(td_div,beta_div,krel,robot_new)
% Performs data set scaling and returns a new vector of velocity magnitudes
%
% Inputs : 
%          - td_div   spacing of td
%          - beta_div spacing of beta
%          - krel         new relative stiffness value
%          - robot_new    robot data structure must have leg length
%
% Outputs:
%          - beta_vec     beta vector
%          - td_vec       touchdown vector
%          - V0           velocity magnitudes
%          - beta_step    returns beta_div
%          - td_step      returns td_div


load("krel_slope_new.mat",'td_vec','beta_vec','Ax_fr','bx_fr');

td_step=td_div;
beta_step=beta_div;
num=abs(ceil((beta_vec(end)-beta_vec(1))/beta_div))+1;
beta_vec_new = linspace(beta_vec(1),beta_vec(end),num);
td_vec_new   = linspace(td_vec(1),td_vec(end),num);


Ax=zeros(length(td_vec_new),length(beta_vec_new));
Bx=zeros(length(td_vec_new),length(beta_vec_new));


for i=1:length(td_vec_new)

    for j=1:length(beta_vec_new)
        
        Ax(i,j)=interp2(beta_vec,td_vec,Ax_fr,beta_vec_new(j),td_vec_new(i));
        Bx(i,j)=interp2(beta_vec,td_vec,bx_fr,beta_vec_new(j),td_vec_new(i));

    end

end

%Calculate matrix of froude values
Fr=Ax*krel+Bx;

%Undo Froudes and convert to magnitudes
V0 = sqrt(Fr*9.81*robot_new.l0);
V0=real(V0);


end