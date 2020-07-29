function [bend,twist,shear,stretch] = four_energies(zvec)
% THIS FILE ONLY FUNCTIONS FOR A SEMICIRCULAR BOUNDARY CONDITIONS
% This file now has self contact gradient terms.
% This file neither constructs nor returns a hessian

global avgs_for_mer
global stiffs_for_mer
global r0
global rn
global q0
global qn
s = size(zvec);
zlen = s(1);
grad = zeros(zlen,1);
nbp = zlen/7+2; %zln/7 = #of unknown frames. total here is 2 knowns + th zlen/7 unknowns
penalty_weight=100;
njunc = nbp-1;


% Setup for building a sparse Hessian (ntriplies = number of times we'll put stuff in the Hessian)
ntriplets = 457*(njunc-2)+163;
I = zeros(ntriplets,1); J = zeros(ntriplets,1); X = zeros(ntriplets,1);

q = zeros(4,nbp);
r = zeros(3,nbp);
id = eye(4);

%assign values to the vectors q and r for all i along the rod
for i=2:nbp-1
    q(:,i)=zvec(4*(i-2)+1:4*(i-1),1);
    r(:,i)=zvec(4*(nbp-2)+3*(i-2)+1:4*(nbp-2)+3*(i-1),1);
end

%impose the boundary conditions.
r(:,end) = rn';
r(:,1) = r0';
q(:,end) = qn';
q(:,1) = q0';


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
        cay(k,i) = 2/dfac(i)*q(:,i+1)'*(bk*q(:,i));
        
    end
    dr = r(:,i+1)-r(:,i);
    dirs = compute_ds(q(:,i+1)+q(:,i));
    tr(1,i) = dr'*dirs(:,1);
    tr(2,i) = dr'*dirs(:,2);
    tr(3,i) = dr'*dirs(:,3);
end


        
%disp(tr)

%energy function and gradients
%also define dj and dj derivatives
bend = 0.0; twist = 0.0; shear = 0.0; stretch = 0.0;
for i=1:njunc
    bend = bend+stiffs_for_mer(i,4)*(cay(1,i)-avgs_for_mer(i,4))^2/2;
    bend = bend+stiffs_for_mer(i,5)*(cay(2,i)-avgs_for_mer(i,5))^2/2;
    twist = twist+stiffs_for_mer(i,6)*(cay(3,i)-avgs_for_mer(i,6))^2/2;
    shear = shear+stiffs_for_mer(i,1)*(tr(1,i)-avgs_for_mer(i,1))^2/2;
    shear = shear+stiffs_for_mer(i,2)*(tr(2,i)-avgs_for_mer(i,2))^2/2;
    stretch = stretch+stiffs_for_mer(i,3)*(tr(3,i)-avgs_for_mer(i,3))^2/2;
end

