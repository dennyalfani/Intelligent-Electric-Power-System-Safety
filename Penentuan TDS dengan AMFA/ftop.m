function top = ftop(x) % x is TMS for each relay
relayIND = [1 2; % relay ID. ncases x 2
            2 0];
Ifault = [13240 1760; % fault felt by relayIND. ncases x 2
          38490    0];
Ipickup = [1820 243]; % I pickup
% CTratio = [1200/5 500/5]; % CT ratio
% V       = [20 150]; % relay voltage level
% Vbase   = 150; %base V
a = 13.5/1.5; % coefficient
b = 1; % expo
%% fault is primary an pickup in primary, so no change
ncases = size(Ifault,1);
% Iused = zeros(ncases,2);
% for i = 1:ncases
%     for j = 1:2
%         relayName = relayIND(i,j);
%         if relayName ~= 0
%             Iused(i,j) = Ifault(i,j)/CTratio(relayName);
%         end
%     end
% end
Iused = Ifault;
%% operating time at each relay each role
waktu = zeros(ncases,2);
for i = 1:ncases
    for j = 1:2
        relayName = relayIND(i,j);
        if relayName ~= 0
            rasioI = Iused(i,j)/Ipickup(relayName);
%             if rasioI > 20 %% comment line 31 - 33 if relay has no max. tap
%                 rasioI = 20;
%             end
            dummy = ((rasioI)^b)-1;
            waktu(i,j) = a*x(relayName)/dummy;
        end
    end
end
disp(waktu);
%% constraint
langgar = 0;
% CTI
CTI = 0.2;
for i = 1:ncases
    if waktu(i,2) - waktu(i,1) < CTI && relayIND(i,2) ~= 0
        langgar = langgar+2000;
    end
end
% operating time
for i = 1:ncases
    for j = 1:2
        relayName = relayIND(i,j);
        if relayName ~= 0
            if waktu(i,j) < 0.1 || waktu(i,j)> 2.5
                langgar = langgar+2000;
            end
        end
    end
end
%% final value
top = sum(sum(waktu))+langgar;
end