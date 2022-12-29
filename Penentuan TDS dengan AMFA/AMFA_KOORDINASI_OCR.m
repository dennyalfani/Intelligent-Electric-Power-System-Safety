close all
clc
clear

format long

%Input
jumlah_relay = input ('masukkan jumlah relay           :');
kv = input('masukkan tegangan relay         :');
kvbase = input('masukkan kv base                :');
FLAmasukan = input('masukkan fla                    :');
Isc_max_prim = input('masukkan isc max primer         :');
Isc_max_back = input('masukkan isc max backup         :');
nct = input('masukkan primer CT              :');
TDSmin = 0.025 ;
TDSmax = 1.2 ;
TDSStep= 0.025;
Target_CTI = input('masukkan Target CTI             :');
Topmin = input('masukkan waktu target minimum   :');
%Topmax = input('masukkan waktu target maksimum  :');


relay_backup_primer = input('masukkan pasangan relay         :');
pasangan_relay = size(relay_backup_primer);%brs kolom
Pair_relay = pasangan_relay(1); %diambil barisnya aja



for i = 1:jumlah_relay
    FLA(i) = ((FLAmasukan(1,i)*kv(1,i))/kvbase);

end
for i = 1:jumlah_relay
    Ipickup(i) = 1.05*FLA(1,i);

end

for i = 1:jumlah_relay
    Isc_max_primer(1,i) = ((Isc_max_prim(1,i)*kv(1,i))/kvbase);
end

for i = 2:jumlah_relay

    Isc_max_backup(1,1) = inf;

    Isc_max_backup(1,i) = ((Isc_max_back(1,i)*kv(1,i))/kvbase);
end


relaysize = [1 jumlah_relay];
terbaik=100;
part=1;

%% FA Parameter
nVar=length(FLA);                 % Number of Decision Variables = Total rele

VarSize=[1 nVar];                      % Decision Variables Matrix Size

MaxIt=50;         % Maximum Number of Iterations

nPop=30;            % Number of Fireflies (Swarm Size)

gamma=1;            % Light Absorption Coefficient

beta0=2;            % Attraction Coefficient Base Value

alpha=0.9;          % Mutation Coefficient

alpha_damp=0.98;    % Mutation Coefficient Damping Ratio

delta=0.05*(TDSmax-TDSmin);     % Uniform Mutation Range

m=2; %untuk pangkat di variable beta


if isscalar(TDSmin) && isscalar(TDSmax)
    dmax = (TDSmax-TDSmin)*sqrt(nVar);
else
    dmax = norm(TDSmax-TDSmin);
end


%Firefly Template
firefly.TDS=[];
firefly.TimeOperationPrim=[];
firefly.TimeOperationSek=[];
firefly.Cost=[];
firefly.CTI.Cost = [];
firefly.CTI.Error = [];

%Global Best
Best_Solution.Cost = inf;

%Create Population
for i=1:nPop
    firefly(i).TDS=ceil(unifrnd(TDSmin,TDSmax,VarSize)/TDSStep)*TDSStep; % bangkit random dari varmin sampai varmax , baris var size kolom var size

end

for i= 1:nPop
    for j = 1:jumlah_relay
        TDSawal(i,j) = firefly(i).TDS(j);
    end
end

%Time Operation
% Evaluasi Top Relay 1 hingga ke n-1
for i=1:nPop
    for j=1:Pair_relay
        if Isc_max_primer(j)>30*Ipickup(j) && Isc_max_backup(j+1)>30*Ipickup(j+1)
            for k=j:jumlah_relay
                firefly(i).TimeOperationPrim(k)=0.466*firefly(i).TDS(k);

            end

            for l=j:Pair_relay
                firefly(i).TimeOperationSek(l+1)=firefly(i).TimeOperationPrim(l+1);

            end

        else if Isc_max_primer(j)>30*Ipickup(j) && Isc_max_backup(j+1)<30*Ipickup(j+1)
                for k=j:jumlah_relay
                    firefly(i).TimeOperationPrim(k)=0.466*firefly(i).TDS(k);
                    firefly(i).TimeOperationSek(k)=firefly(i).TDS(k)*(13.5/((Isc_max_backup(1,k)/Ipickup(1,k))-1));

                end

        else if Isc_max_primer(j)<30*Ipickup(j) && Isc_max_backup(j+1)<30*Ipickup(j+1)
                for k=j:jumlah_relay
                    firefly(i).TimeOperationPrim(k)=firefly(i).TDS(k)*(13.5/((Isc_max_primer(1,k)/Ipickup(1,k))-1));
                    firefly(i).TimeOperationSek(k)=firefly(i).TDS(k)*(13.5/((Isc_max_backup(1,k)/Ipickup(1,k))-1));

                end
        end
        end
        end
    end

    % Evaluasi Top Relay ke n
    if j==Pair_relay
        for t=j+1
            if Isc_max_primer(t)>30*Ipickup(t)
                firefly(i).TimeOperationPrim(t)=0.466*firefly(i).TDS(t);
                num2str(firefly(i).TimeOperationPrim(t),'%7.10f');
            else if Isc_max_primer(t)<30*Ipickup(t)
                    firefly(i).TimeOperationPrim(t)=firefly(i).TDS(t)*(13.5/((Isc_max_primer(1,t)/Ipickup(1,t))-1));

            end
            end
        end
    end
end

%CTI
for i=1:nPop
    for j=1:jumlah_relay-1
        firefly(i).CTI.Cost(j) = firefly(i).TimeOperationSek(j+1)- firefly(i).TimeOperationPrim(j);
        firefly(i).CTI.Error(j)= firefly(i).CTI.Cost(j) - Target_CTI;

    end
end

%objective function
for i=1:nPop
    firefly(i).Cost = sum(firefly(i).TimeOperationPrim);
end


% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

fireflybaru.TDS = [];
fireflybaru.Cost = [];
fireflybaru.TimeOperationPrim = [];
fireflybaru.TimeOperationSek = [];
fireflybaru.CTI = [];
datafirefly.TDS = [];
datafirefly.TimeOperationPrim = [];
datafirefly.TimeOperationSek = [];

newpop = repmat(firefly,1,1);

%% Firefly Algorithm Main Loop
for it=1:MaxIt

    for i=1:nPop
        for j=1:nPop
            TDSmin = 0.025;
            TDSmax = 1.2;
            if newpop(j).Cost < newpop(i).Cost
                rij=norm(newpop(i).TDS-newpop(j).TDS)/dmax;
                beta=beta0*exp(-gamma*rij^m);
                e=delta*unifrnd(-1,+1,VarSize);


                datafirefly(j,i).TDS = ceil((newpop(i).TDS ...
                    + beta*rand(VarSize).*(newpop(j).TDS-newpop(i).TDS) ...
                    + alpha*e)/TDSStep)*TDSStep;


            else datafirefly(j,i).TDS = firefly(i).TDS;

            end

            %hitung t operation
            for c=1:Pair_relay
                if Isc_max_primer(c)>30*Ipickup(c) && Isc_max_backup(c+1)>30*Ipickup(c+1)
                    for v=c:jumlah_relay
                        datafirefly(j,i).TimeOperationPrim(v)=0.466*datafirefly(j,i).TDS(v);

                    end

                    for l=c:Pair_relay
                        datafirefly(j,i).TimeOperationSek(l+1)=datafirefly(j,i).TimeOperationPrim(l+1);

                    end

                else if Isc_max_primer(c)>30*Ipickup(c) && Isc_max_backup(c+1)<30*Ipickup(c+1)
                        for v=c:jumlah_relay
                            datafirefly(j,i).TimeOperationPrim(v)=0.466*datafirefly(j,i).TDS(v);
                            datafirefly(j,i).TimeOperationSek(v)=datafirefly(j,i).TDS(v)*(13.5/((Isc_max_backup(1,v)/Ipickup(1,v))-1));

                        end

                else if Isc_max_primer(c)<30*Ipickup(c) && Isc_max_backup(c+1)<30*Ipickup(c+1)
                        for v=c:jumlah_relay
                            datafirefly(j,i).TimeOperationPrim(v)=datafirefly(j,i).TDS(v)*(13.5/((Isc_max_primer(1,v)/Ipickup(1,v))-1));
                            datafirefly(j,i).TimeOperationSek(v)=datafirefly(j,i).TDS(v)*(13.5/((Isc_max_backup(1,v)/Ipickup(1,v))-1));

                        end
                end
                end
                end
            end

            % Evaluasi Top Relay ke n
            if c==Pair_relay
                for t=c+1
                    if Isc_max_primer(t)>30*Ipickup(t)
                        datafirefly(j,i).TimeOperationPrim(t)=0.466*datafirefly(j,i).TDS(t);

                    else if Isc_max_primer(t)<30*Ipickup(t)
                            datafirefly(j,i).TimeOperationPrim(t)=datafirefly(j,i).TDS(t)*(13.5/((Isc_max_primer(1,t)/Ipickup(1,t))-1));

                    end
                    end
                end
            end



            %Error CTI
            for k=1:jumlah_relay-1

                datafirefly(j,i).CTI.Cost(k) = datafirefly(j,i).TimeOperationSek(k+1)- datafirefly(j,i).TimeOperationPrim(k);
                datafirefly(j,i).CTI.Error(k)= datafirefly(j,i).CTI.Cost(k) - Target_CTI;

            end

            %Objective function evaluation
            datafirefly(j,i).Cost = sum(datafirefly(j,i).TimeOperationPrim);



            %Error Cost Evaluation
            for k=1:jumlah_relay
                if datafirefly(j,i).TimeOperationPrim(k)<Topmin %|| datafirefly(j,i).TimeOperationPrim(k)>Topmax
                    datafirefly(j,i).Cost=100;
                end
            end

            for k=1:Pair_relay
                if datafirefly(j,i).CTI.Cost(k)<Target_CTI
                    datafirefly(j,i).Cost=100;
                end
            end




        end

    end

    %update personal best
    Best.Cost = inf;
    for y = 1:nPop
        Best.Cost = inf;
        for u = 1:nPop
            if datafirefly(u,y).Cost<Best.Cost;
                newpop(y) = datafirefly(u,y);

                Best = datafirefly(u,y);
                %
            end
        end
    end


    for j=1:nPop
        datanewpop(j,it)= newpop(j);%data newpop tiap iterasi
    end

    %update global best
    for y = 1:nPop
        if newpop(y).Cost < Best_Solution.Cost
            Best_Solution  = newpop(y);
        end
    end

    for i=1:jumlah_relay
        BestSolit(it,i)=Best_Solution.TDS(i);
    end


    %besttds and cost per iterasi
    for z=1:jumlah_relay
        BestTDS(it,z) = Best_Solution.TDS(z);
    end
    BestCost(it)=Best_Solution.Cost;


    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);

    % Damp Mutation Coefficient

    alpha = alpha*((1/(2*MaxIt))^(1/(MaxIt+1)));
    %alpha = alpha - (0.02*(it));
    %alpha = 0.9;
end


for i=1:MaxIt-1
    for j=1:jumlah_relay
        Deltaerror(i,j) = BestSolit(i+1,j)-BestSolit(i,j);
    end
end

for i = 1:jumlah_relay
    pickup(i) = (ceil((Ipickup(1,i)/nct(1,i))/0.01)*0.01)*kvbase/kv(1,i);

end

for i = 1:jumlah_relay
    nilaipickup(i) = pickup(1,i)*nct(1,i);


end

disp(' Relay No. | pickup | Ipickup');

for m = 1:jumlah_relay
    fprintf('  %5.3g', m);
    fprintf('  %8.2f', pickup(m));
    fprintf('  %10.2f', nilaipickup(m));
    fprintf('\n');


end

disp('                            TIME DIAL                             ');
disp('==================================================================');
disp('          | Relay No. |   TDS  |  WAKTU OPERASI |  WAKTU OPERASI |');
disp('                               |     PRIMER     |    SEKUNDER    |');
disp('------------------------------------------------------------------');
for m = 1:jumlah_relay
    fprintf('  %13.3g', m);
    fprintf('  %13.4f', Best_Solution.TDS(m));
    fprintf('  %10.4f', Best_Solution.TimeOperationPrim(m));
    fprintf('  %14.4f', Best_Solution.TimeOperationSek(m));
    fprintf('\n');
end


figure
title('Best Cost')
plot(BestCost, 'LineWidth', 2);
xlabel('Iterations');
ylabel('Best Cost');
grid on

for j=1:jumlah_relay
    figure;
    plot(BestTDS(:,j),'LineWidth',2);
    title(['TDS Rele('  num2str(j)   ') ']);
end

for j = 1:jumlah_relay
    figure;
    plot(TDSawal(:,j),'o');
    title(['Persebaran Awal TDS Rele('  num2str(j)   ') ']);
end

for j = 1:jumlah_relay
    figure;
    plot(BestSolit(:,j),'o');
    title(['Pergerakan TDS Rele('  num2str(j)   ') ']);
end

for j = 1:jumlah_relay
    figure;
    plot(abs(Deltaerror(:,j)));
    title(['Delta Error TDS Rele('  num2str(j)   ') '])
end

for j = 1:jumlah_relay
    figure;
    plot((abs(Deltaerror(:,j))),'o');
    title(['Delta Error TDS Rele('  num2str(j)   ') '])
end




%PLOT TCC
figure;

iscx1=[Ipickup(1):1:(30*Ipickup(1))]; %Elemen Sumbu X
iscx2=[Ipickup(2):1:(30*Ipickup(2))]; %Elemen Sumbu X
y1=(13.5*(BestTDS(50,1)))./((iscx1./(Ipickup(1)))-1); %Elemen Sumbu Y
y2=(13.5*(BestTDS(50,2)))./((iscx2./(Ipickup(2)))-1); %Elemen Sumbu Y

kurvatcc=loglog(iscx1,y1,iscx2,y2); %plot TCC
set(kurvatcc,'linewidth',1);
hold on;

kurvaiscprim1=loglog([Isc_max_primer(1) Isc_max_primer(1)],[.01 1e3],'--r'); %PlotIsc Primer
set(kurvaiscprim1,'linewidth',1.5);
hold on;

kurvaiscsek1=loglog([Isc_max_backup(2) Isc_max_backup(2)],[.01 1e3],'--b'); %Plot Isc Backup
set(kurvaiscsek1,'linewidth',1.5);
hold on;

kurvasat=loglog([30*Ipickup(1) (max(Isc_max_primer(2)+20000))],[((BestTDS(50,1))*0.466) ((BestTDS(50,1))*0.466)],[30*Ipickup(2) (max(Isc_max_primer(2)+20000))],[((BestTDS(50,2))*0.466) ((BestTDS(50,2))*0.466)]); %Plot Saturasi
set(kurvasat,'linewidth',1);
hold on
grid on

legend('Kurva Relay Primer','Kurva Relay Backup','Isc Primer 1','Isc Backup 1');
set(kurvatcc(1),'Color','red');
set(kurvasat(1),'Color','red');
set(kurvatcc(2),'Color','blue');
set(kurvasat(2),'Color','blue');

xlim([0,max((Isc_max_primer(2)+20000))]);
ylim([0,1e3]);
title('Kurva TCC');
xlabel('Arus (A)');
ylabel('Waktu (s)');

ax=gca;
ax.XTick=[0.01 300 500 1000 3000 5000 10000 30000 50000 100000];
ax.YTick=[0.1 0.3 0.5 1 3 5 10 30 50 100 300 500 1000];