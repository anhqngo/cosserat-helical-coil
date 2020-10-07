function energy = helix_energy(nseg, stretch_param, contact_strength,...
                                R, c, tmax)
    % The formula for the helix is r(t) = [R*cos(t), R*sin(t), c*t]
    % The following block of code sets the default parameters
    arguments
        nseg (1,1) double = 64;
        stretch_param (1,1) double = 0 % the range is (-1,1)
        contact_strength (1,1) double = 0
        R (1,1) double = 1
        c (1,1) double = 0.1
        tmax (1,1) double = 6*pi
    end
    
    PRINT=false; % parameter to print out the helix
    
    global avgs_for_mer
    global stiffs_for_mer
    global Q K
    global r0
    global rn
    global q0
    global qn
    
    % Set up variables for easier reference
    denom = sqrt(R^2+c^2);
    arcLength = denom*tmax+1e-7; % perturb to avoid MATLAB fminunc fails
    dt = arcLength/nseg;
    total_height = c*tmax;
    displacement = total_height*stretch_param;    
    
    % Parameters for the energy
    contact_range_param = 1;
    Q = contact_strength/nseg^2;
    K = contact_range_param/2;
    
    avgs_for_mer = zeros(nseg+1,6);
    avgs_for_mer(:,3) = (arcLength)/nseg;
    avgs_for_mer(:,4) = 0;
    avgs_for_mer(:,5) = R/(R^2+c^2)*(arcLength/nseg);
    avgs_for_mer(:,6) = c/(R^2+c^2)*(arcLength/nseg);    
    
    stiffs_for_mer = zeros(nseg+1,6);
    k1 = 1;
    k2 = k1;
    k3 = 100;
    a1 = 300;
    a2 = a1;
    a3 = 3000; 
    for i=1:nseg+1
        stiffs_for_mer(i,1:3) = [a1 a2 a3]*nseg;
        stiffs_for_mer(i,4:6) = [k1 k2 k3]*nseg;
    end
            
    % Inextensible formula - compute qs and rs
    rs = zeros(nseg+1,3);
    qs = zeros(nseg+1,4);
    for i = 1:nseg+1
       length = (i-1)*dt;
       s1 = sin(length/denom);   
       c1 = cos(length/denom);
       minus_var = sqrt((denom-c)*(c1+1)/denom);
       plus_var = sqrt((denom+c)*(c1+1)/denom);

       rs(i,:) = [R*c1, R*s1, c*length/denom];
       qs(i,:) = [-(R*s1)/(2*denom*plus_var),...
                  1/2*minus_var,...
                  1/2*plus_var,...
                  -(R*s1)/(2*denom*minus_var)];
    end
    r0 = rs(1,:);
    rn = rs(nseg+1,:);
    rn(3) = rn(3)+displacement/2;
    r0(3) = r0(3)-displacement/2;
    q0 = qs(1,:);
    qn = qs(nseg+1,:);
    
    % Assemble all but the first and last r's and q's into zvecs (vector of
    % variables)
    zvecs = zeros(7*(nseg-1),1);
    rTemp = rs(2:end-1,:);
    qTemp = qs(2:end-1,:);
    for i=1:nseg-1
        zvecs(4*(nseg-1)+3*(i-1)+1:4*(nseg-1)+3*i,1) = rTemp(i,:)';
        zvecs(4*(i-1)+1:4*i,1) = qTemp(i,:)';
    end
    
    % Minimize the discrete-rod energy
    options = optimset('Display','off','MaxIter',4000,'GradObj','on',...
                        'Hessian','off','TolFun',1e-6, 'HessUpdate',...
                         'bfgs','MaxFunEvals',1e+5);
    [zeq,~,~,~] = fminunc('discrete_dna_penalty_en_grad_s',zvecs,options);

    energy = discrete_dna_penalty_en_grad_s(zeq);

    if PRINT
        fprintf("The energy is %f\n", energy)
        disp("Four components (4th scaled by 1000) of energy are");
        [e1,e2,e3,e4]=helix_four_energies(zeq);
        [e1 e2 e3 1000*e4]

        % Save the solution to a file
        filename = 'you_choose.txt';
        fileID2 = fopen(filename,'wt');
        fprintf(fileID2,'%f\r\n',zeq(:,:));
        fclose(fileID2);

        % ---------------------------------------------------------------------

        newNSeg = nseg-1;
        newQs = zeros(newNSeg, 4);
        newRs = zeros(newNSeg, 3);
        fileID = fopen(filename,'r');
        formatSpec = '%f';
        A = fscanf(fileID,formatSpec);
        entry = 1;
        for i = 1:newNSeg
            for j = 1:4
                newQs(i+1,j) = A(entry);
                entry = entry + 1;
            end
        end
        for i = 1:newNSeg
            for j = 1:3
                newRs(i+1,j) = A(entry);
                entry = entry + 1;
            end
        end
        newRs(1,:) = r0;
        newRs(newNSeg+2,:) = rn;
        plot3(newRs(:,1),newRs(:,2),newRs(:,3),'LineWidth',2);
        axis equal;
        grid on;
    end
end