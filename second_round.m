%--------------------------------------------------------------------------
% TOTAL ENERGY WITH DIFFERENT DISPLACEMENT PARAMS + CONTACT STRENGTH PARAMS
%--------------------------------------------------------------------------
plot=true;

if plot
NSEGS = 256;
STRETCH_PARAMS = -[0.03333 0.06667 0.1 0.16667 0.20 0.23333 0.3 0.33333 0.36667];
contact_strength_params = [0, 800, 2000, 3000];
stretch_energies = zeros(length(contact_strength_params), length(STRETCH_PARAMS));

figure(1)
for j=1:length(contact_strength_params)
    for i=1:length(STRETCH_PARAMS)
        stretch_energies(j,i) = helix_energy(NSEGS,STRETCH_PARAMS(i), contact_strength_params(j));
        temp = [STRETCH_PARAMS(i), contact_strength_params(j), stretch_energies(j,i)];
        fprintf('%s\n', sprintf('%d ', temp));       
    end
%     txt = sprintf('Contact strength=%0.2f', contact_strength_params(j)/(64^2));
%     plot(STRETCH_PARAMS, stretch_energies(j,:), 'DisplayName', txt);
%     hold on;
end
% hold off;
% xlabel('Stretch parameters');
% xticks(STRETCH_PARAMS);
% ylabel('Total energy');
% legend('Location','southwest');
% title('Energies of stretched rods with different stretching parameters');
end

%----------------------------------------------------------
% CONVERGENCE OF TOTAL ENERGY WITH DIFFERENT NSEG
%----------------------------------------------------------
plot=false;
if plot
NSEGS = [16, 32, 64, 128, 256];
ENERGIES = [0, 0, 0, 0, 0];

for j = 1:length(NSEGS)
    ENERGIES(j) = helix_energy(NSEGS(j));
end
semilogx(NSEGS, ENERGIES);
xticks(NSEGS);
ylabel("Total nergy");
xlabel('Number of segments');
end


%----------------------------------------------------------
% CONVERGENCE OF ENERGY WITH DIFFERENT NSEG
%----------------------------------------------------------
plot = false;
if plot
NSEGS = [16, 32, 64, 128, 256];
ENERGIES = [27.8085    0.0426    0.0313    3.0706;...
            9.3984    0.0152    0.0122    1.1710; ...
            2.4929    0.0220    0.0105    0.6999; ...
            0.6039    0.0053    0.0025    0.1669; ...
            0.2509    0.0010    0.0024    0.3767];
ENERGY_NAMES = ["Bending","Twisting","Shearing","Stretching"];

figure(1)
for j = 1:4
%     plot(NSEGS, ENERGIES(:,j), 'DisplayName',ENERGY_NAMES(j));
    loglog(NSEGS, ENERGIES(:,j), '-s','DisplayName',ENERGY_NAMES(j));
    hold on;
end
xticks(NSEGS)
ylabel("Energy");
legend('Location','northeast');
xlabel('Number of segments');
end


% 
% % ----------------------------------------------------------------------
% figure(2)
% for j=1:m
%     txt = sprintf('Contact strength=%0.2f', contact_strength_params(j)/(64^2));
%     plot(stretch_params, stretch_energies(j,:), 'DisplayName', txt);
%     hold on;
% end
% hold off;
% xlabel('Stretch parameters');
% ylabel('Total energy');
% legend('Location','northeast');
% title('Total stretch energy with different stretching parameters');
% 
% 
% % ----------------------------------------------------------------------
% figure(3)
% for j=1:m
%     txt = sprintf('Contact strength=%0.2f', contact_strength_params(j)/(64^2));
%     plot(stretch_params, shrink_energies(j,:), 'DisplayName', txt);
%     hold on;
% end
% hold off;
% xlabel('Shrink parameters');
% ylabel('Total energy');
% legend('Location','northwest');
% title('Total shrink energy with different shrinking parameters');
% 
