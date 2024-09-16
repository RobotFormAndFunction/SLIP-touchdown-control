function [magsol]=feasible_points(td_vec,beta_vec,magsol,robot)
% Removes all initial conditions that won't allow system to clear leg
% length while jumping
%
% Inputs : 
%          - td_vec   td vector
%          - beta_vec beta vector
%          - magsol   initial velocity matrices
%          - robot    robot data structure must have leg length
%
% Outputs:
%          - V0           velocity magnitudes


        %Check if system can clear step
c=0;        

for i=1:length(td_vec)
        TO_angle=td_vec(i);
    for j=1:length(beta_vec)

        if magsol(i,j)~=0

            yvel_TO=magsol(i,j)*sin(beta_vec(j));

            if (yvel_TO^2-2*9.81*(robot.l0-robot.l0*sin(TO_angle)))<0
                magsol(i,j)=0;
                c=c+1;
            end
        end
    end
end


end