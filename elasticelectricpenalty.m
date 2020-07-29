function [elastic1,elastic2,electric,penalty] = elasticelectricpenalty(zvec)
% THIS FILE ONLY FUNCTIONS FOR A SEMICIRCULAR BOUNDARY CONDITIONS
% This file now has self contact gradient terms.
% This file neither constructs nor returns a hessian

global avgs_for_mer
global stiffs_for_mer
s = size(zvec);
zlen = s(1);

nbp = zlen/7+2; %zln/7 = #of unknown frames. total here is 2 knowns + th zlen/7 unknowns
penalty_weight=100;
njunc = nbp-1;

q = zeros(4,nbp);
r = zeros(3,nbp);

%assign values to the vectors q and r for all i along the rod
for i=2:nbp-1
    q(:,i)=zvec(4*(i-2)+1:4*(i-1),1);
    r(:,i)=zvec(4*(nbp-2)+3*(i-2)+1:4*(nbp-2)+3*(i-1),1);
end

%institute the semi-circular boundary conditions.
global endpoint twistAngle 
r(:,end) = [3 ; 0 ; 0];
r(:,1) =  [3 ; 0 ; 0];
%the Cayley form %make sure this is -1 for the x-Z semi-circle.
q(:,1) = [-0.5257 ; 0 ; 0 ; 0.8507];
q(:,nbp) = [0.5257 ; 0 ; 0 ; -0.8507];
%totalturns = sum(avgs_for_mer(1:nbp-1,6)/5)/(2*pi);
% totalturns_in_ddpegh = totalturns
% intturns_in_ddpegh = intturns
%either same as beginning quaternion or negative of it.

%twisted rod
%needs to be generalized for twisting greater than 2pi
%if twistAngle >=0 && twistAngle <= 1*pi
%    q(:,nbp) = [sqrt((1-cos(twistAngle))/2),sqrt((1+cos(twistAngle))/2),0,0];
%elseif twistAngle >= pi && twistAngle <= 2*pi
%    q(:,nbp) = [sqrt((1-cos(twistAngle))/2),-sqrt((1+cos(twistAngle))/2),0,0];
%elseif twistAngle >= 2*pi && twistAngle <= 3*pi
%     q(:,nbp) = [-sqrt((1-cos(twistAngle))/2),-sqrt((1+cos(twistAngle))/2),0,0];
%elseif twistAngle >= 3*pi && twistAngle <= 4*pi
%    q(:,nbp) = [-sqrt((1-cos(twistAngle))/2),sqrt((1+cos(twistAngle))/2),0,0];
%elseif twistAngle >= 4*pi && twistAngle <= 5*pi
%    q(:,nbp) = [sqrt((1-cos(twistAngle))/2),sqrt((1+cos(twistAngle))/2),0,0];
%elseif twistAngle >= 5*pi && twistAngle <= 6*pi
%    q(:,nbp) = [sqrt((1-cos(twistAngle))/2),-sqrt((1+cos(twistAngle))/2),0,0];
%elseif twistAngle >= 6*pi && twistAngle <= 7*pi
%     q(:,nbp) = [-sqrt((1-cos(twistAngle))/2),-sqrt((1+cos(twistAngle))/2),0,0];
%elseif twistAngle >= 7*pi && twistAngle <= 8*pi
%    q(:,nbp) = [-sqrt((1-cos(twistAngle))/2),sqrt((1+cos(twistAngle))/2),0,0];
%end


%define the B matrices used in theta and a definitions
b = zeros(3,4,4);
b(1,1,4) = 1; b(1,2,3) = 1; b(1,3,2)=-1; b(1,4,1) = -1;
b(2,1,3) = -1; b(2,2,4) = 1; b(2,3,1) = 1; b(2,4,2) = -1;
b(3,1,2) = 1; b(3,2,1) = -1; b(3,3,4) = 1; b(3,4,3) = -1;
bq = zeros(3,4,nbp);
for i=1:3
    bi = zeros(4,4);
    for j1 = 1:4
        for j2 = 1:4
            bi(j1,j2)=b(i,j1,j2);
        end
    end
    for k = 1:nbp
        dum = bi*q(:,k);
        for j1=1:4
            bq(i,j1,k)=dum(j1);
        end
    end
end

%cay = thetas
%tr = a's (shift/slide/rise)
%dfac = qa * qb (q dot its neighbor)

cay = zeros(3,njunc);
tr = zeros(3,njunc);
dfac = zeros(1,njunc);
for i=1:njunc
    dfac(i) = q(:,i+1)'*q(:,i);
    for k=1:3
        bk = zeros(4,4);
        for j1 = 1:4
            for j2 = 1:4
                bk(j1,j2)=b(k,j1,j2);
            end
        end
        cay(k,i) = 10/dfac(i)*q(:,i+1)'*(bk*q(:,i));
        
    end
    dr = r(:,i+1)-r(:,i);
    dirs = compute_ds(q(:,i+1)+q(:,i));
    tr(1,i) = dr'*dirs(:,1);
    tr(2,i) = dr'*dirs(:,2);
    tr(3,i) = dr'*dirs(:,3);
end


%energy function and gradients
%also define dj and dj derivatives
elastic1 = 0.0;
elastic2 = 0.0;
penalty = 0.0;
for i=1:njunc
    for j1 = 1:3
        elastic1 = elastic1+stiffs_for_mer(i,j1+3)*(cay(j1,i)-avgs_for_mer(i,j1+3))^2/2;
        elastic2 = elastic2+stiffs_for_mer(i,j1)*(tr(j1,i)-avgs_for_mer(i,j1))^2/2;
    end
    penalty = penalty+penalty_weight*(q(:,i+1)'*q(:,i+1)-1)^2;
end

% self contact energy portion
global Q K
electric = 0;
xnorm = 0.0;
for i=1:nbp
    for j=i:nbp
%        if i==1 && j==nbp
            %nothing
%        elseif abs(i-j) > nbp/10
        if abs(i-j) > nbp/10 && abs(i+(nbp-j)) > nbp/10
            xnm = r(:,i) - r(:,j);
            xnorm = norm(xnm,2);
            electric = electric + exp(-K*xnorm) / xnorm;
        end
    end
end

electric = Q*electric;
