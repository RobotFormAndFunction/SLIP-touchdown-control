function [x,dx,y,dy,ddx,ddy]=p2r(p, dir)
%p2r converts from polar cordinates and velocities to rectangular
%note theta 0 is 0 on when || to the the x axis in the third quatrant and
%clockwise is the + theta direction
% X,Y are defined as normal cartesian directions
%Inputs [p], polar vector with fields:
%          [r,dr/dt,theta,dtheta/dt]
%          if given a nx6 matrix will return accelerations in that case it
%          should be formated as:
%          [r,r/dt,r/ddt,theta,theta/dt,theta/ddt]
%
%        dir, direction optional. Enter a 1 if the input is swapped
%             [theta,theta/dt,theta/ddt,r,r/dt,r/ddt]
%
%Outputs [x,dx/dt,y,dy/dt, dx/ddt,dy/ddt]

if ~exist('dir','var') dir=0;  end

x=zeros(length(p(:,1)),1);
y=x;
dx=x;
dy=y;
ddx=x;
ddy=y;
c=p;


if length(p(1,:))==4
    for i=1:length(p(:,1))
        if dir==1
          p=[c(i,3) c(i,4) c(i,1) c(i,2)];  
        else
          p=c(i,:);  
        end
        
        x(i)=-p(1)*cos(p(3));
        y(i)=p(1)*sin(p(3));
        
        dx(i)=-p(2)*cos(p(3))+p(1)*p(4)*sin(p(3));
        dy(i)=p(2)*sin(p(3))+p(1)*p(4)*cos(p(3));
    end
end

if length(p(1,:))==6
    for i=1:length(p(:,1))
        
        if dir==1
          p=[c(i,4) c(i,5) c(i,6) c(i,1) c(i,2) c(i,3)];  
        else
          p=c(i,:);  
        end
        
        x(i)=-p(1)*cos(p(4));
        y(i)=p(1)*sin(p(4));
        
        dx(i)=-p(2)*cos(p(4))+p(1)*p(5)*sin(p(4));
        dy(i)=p(2)*sin(p(4))+p(1)*p(5)*cos(p(4));
        
        ddx(i)=(p(1)*p(5)^2-p(3))*cos(p(4))+(p(1)*p(6)+2*p(2)*p(5))*sin(p(4));
        ddy(i)=(p(3)-p(1)*p(5)^2)*sin(p(4))+(p(1)*p(6)+2*p(2)*p(5))*cos(p(4));
    end
end



end