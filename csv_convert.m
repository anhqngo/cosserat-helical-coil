M = readmatrix('matlab_runs.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
for k=0:2
    start_idx = 1+13*k+52;
    end_idx = 13*(k+1)+52;
    stretch_params = M(start_idx:end_idx, 1);
    energy = M(start_idx:end_idx, 3);
    contact_strength_param = M(start_idx, 2);
    [stretch_params, sorted_order] = sort(stretch_params, 'ascend');
    energy = energy(sorted_order);    
    
    txt = sprintf('Contact strength=%0.2f', contact_strength_param/(256^2));    

    p = polyfit(stretch_params, energy, 2);
    
    % find the minimum
    d1p = polyder(p);                           % First Derivative
    d2p = polyder(d1p);                         % Second Derivative
    ips = roots(d1p);                           % Inflection Points
    xtr = polyval(d2p, ips);                    % Evaluate ‘d2p’ at ‘ips’
    x0 = ips((xtr > 0) & (imag(xtr)==0));       % Find Minima
    y0 = polyval(p, x0);
    
    L = 18.9436;
    
    k = (energy - y0)*2./(((stretch_params-x0)*L).^2);
    hold on;
    txt = sprintf('Contact strength=%0.2f', contact_strength_param/(256^2));    
    scatter(stretch_params, k, 'DisplayName', txt);
    
end
hold off;
xlabel("Stretching parameters");
ylabel("k value");
ylim([-2,2]);
    set ( gca, 'xdir', 'reverse' );

legend('Location','southwest');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if false
figure(1);
for k=0:3
    start_idx = 1+13*k;
    end_idx = 13*(k+1);
    stretch_params = M(start_idx:end_idx, 1);
    energy = M(start_idx:end_idx, 3);
    contact_strength_param = M(start_idx, 2);

    txt = sprintf('Contact strength=%0.2f', contact_strength_param/(256^2));    
    
    [stretch_params, sorted_order] = sort(stretch_params, 'descend');
    energy = energy(sorted_order);
    
    disp(stretch_params);
    
    plot(stretch_params, energy, 'DisplayName', txt);
    set ( gca, 'xdir', 'reverse' )
    hold on;
end

hold off;
xlabel('Stretch parameters');
ylabel('Total energy');
legend('Location','southwest');
title('Energies of stretched rods with different stretching parameters');
end