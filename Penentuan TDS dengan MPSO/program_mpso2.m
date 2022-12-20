%% ===============Relays Data=============== %%
% Saving The TDS
for i=1:nPop
    Particle(i).xTDS(1)=GlobalBest.xTDS(2);
    Particle(i).TimeOperationPrim(1)=GlobalBest.TimeOperationPrim(2);
end

for j=1:main_backup_pair
    Target_TopMin=GlobalBest.TimeOperationPrim(j+1);
    Iscmax_prim(1,j)=Iscmax_prim(j+1);
    Iscmax_back(1,j)=Iscmax_back(j+1);
    kV(1,j)=kV(j+1);
    FLA(1,j)=FLA(j+1);
    CT(1,j)=CT(j+1);
    ConstantSat(j)=ConstantSat(j+1);
    CurveType(j)=CurveType(j+1);
    coef_k(j)=coef_k(j+1);
    coef_alpha(j)=coef_alpha(j+1);
    coef_beta(j)=coef_beta(j+1);
end

% Input Relays Data
fprintf('\n');
fprintf('Masukkan Tegangan Relay ke %.0f (kV)',part);
kV(1,2)=input(' : ');
fprintf('Masukkan Primer CT ke %.0f ',part);
CT(1,2)=input(' : ');
fprintf('Masukkan FLA relay ke %.0f (A)',part);
FLA(1,2)=input(' : ');
fprintf('Masukkan Iscmax Relay Primer ke %.0f (A)',part);
Iscmax_prim(1,2)=input(' : ');
fprintf('Masukkkan Iscmax Relay Backup ke %.0f (A)',part);
Iscmax_back(1,2)=input(' : ');
fprintf('Jenis Kurva :\n');
fprintf('(1) Standard Inverse\n');
fprintf('(2) Very Inverse\n');
fprintf('(3) Long Time Inverse\n');
fprintf('(4) Extremely Inverse\n');
fprintf('(5) Ultra Inverse\n');
CurveType(2)=input('Pilih Jenis Kurva : ');

if CurveType(2)==1
    coef_k(2)=0.14;
    coef_alpha(2)=0.02;
    coef_beta(2)=2.97;
else if CurveType(2)==2
        coef_k(2)=13.5;
        coef_alpha(2)=1;
        coef_beta(2)=1.5;
    else if CurveType(2)==3
            coef_k(2)=120;
            coef_alpha(2)=1;
            coef_beta(2)=13.33;
        else if CurveType(2)==4
                coef_k(2)=80;
                coef_alpha(2)=2;
                coef_beta(2)=0.808;
            else if CurveType(2)==5
                    coef_k(2)=315.2;
                    coef_alpha(2)=2.5;
                    coef_beta(2)=1;
                end
            end
        end
    end
end

if kV(2)~=kV(1)
    fprintf('==========================================\n');
    fprintf('                PERHATIAN                  \n');
    fprintf('==========================================\n');
    fprintf('LEVEL TEGANGAN BERBEDA TELAH TERDETEKSI \n');
    fprintf('Apakah CTI Ingin Diminimalkan?\n');
    fprintf('YA = 1, TIDAK = 2\n');
    CTIchoice=input('Masukkan Pilihan : ');
    if CTIchoice==1
        fprintf('STATUS : CTI DIMINIMALKAN\n');
        fprintf('==========================================');
        Target_TopMin=0.1;
        Target_CTI=0;
        %             Target_CTI=0.01;
        fprintf('\n');
    else if CTIchoice==2
            Target_CTI=input('Target CTI (s) : ');
            Target_TopMin=0.1;
            fprintf('\n');
        end
    end
else Target_CTI=input('Target CTI (s) : ');
end

relay_primary_backup=[1 2];

while part<=maxpart
    while CurveType<=5
        % Inisialisasi Ipickup
        for j=1:nRelay
            FLA_relay(1,j)=((FLA(j)*kV(j))/kVbase);
        end
        
        for j=1:nRelay
            Ip(j)=1.05*FLA_relay(j);
        end
        
        for j=1:nRelay
            Iscmax_primer(1,j)=((Iscmax_prim(j)*kV(j))/kVbase);
            Iscmax_backup(1,j)=((Iscmax_back(j)*kV(j))/kVbase);
        end
        
        ipsat=20; %Berapa tap saturasinya
        ConstantSat(2)=((coef_k(2)/coef_beta(2))*1)/((((ipsat/1))^coef_alpha(2))-1);
        
        %% ===============Initialization of PSO===============%%
        GlobalBest.Cost=inf;
        
        % Initialize Population Members
        for i=1:nPop
            % Generate Random Solution of TDS
            Particle(i).xTDS(2)=ceil((((TDSmax-TDSmin)*rand)+TDSmin)*1/stepTDS)/(1/stepTDS);
            
            % Initialize vTDS
            Particle(i).vTDS=0.1*Particle(i).xTDS;
            %Particle(i).vTDS=zeros(RelaySize);
            
        end
        
        % Array to Hold Initial Random Solution of TDS
        for i=1:nPop
            SaveTDS_initial(i,part)=Particle(i).xTDS(2);
        end
        
        for i=1:nPop
            SaveTDS_initial2(1,i+((part-1)*nPop))=Particle(i).xTDS(2);
        end
        
        % Time operation
        for i=1:nPop
            if Iscmax_primer(1)>20*Ip(1) && Iscmax_backup(2)>20*Ip(2)
                Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
                Particle(i).TimeOperationSek(2)=Particle(i).TimeOperationPrim(2);
            else if Iscmax_primer(1)>20*Ip(1) && Iscmax_backup(2)<20*Ip(2)
                    Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
                    Particle(i).TimeOperationSek(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_backup(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                else if Iscmax_primer(1)<20*Ip(1) && Iscmax_backup(2)<20*Ip(2)
                        Particle(i).TimeOperationPrim(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_primer(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                        Particle(i).TimeOperationSek(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_backup(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                    end
                end
            end
            if Iscmax_primer(2)>20*Ip(2)
                Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
            else if Iscmax_primer(2)<20*Ip(2)
                    Particle(i).TimeOperationPrim(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_primer(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                end
            end
        end
        
        for i=1:nPop
            SaveTopPrim_initial(i,part)=Particle(i).TimeOperationPrim(2);
            SaveTopBack_initial(i,part)=Particle(i).TimeOperationSek(2);
        end
        
        
        for i=1:nPop
            while Particle(i).TimeOperationPrim(2)<Target_TopMin
                if TDSmin<=TDSmax
                    TDSmin=TDSmin+stepTDS;
                    Particle(i).xTDS(2)=ceil((((TDSmax-TDSmin)*rand)+TDSmin)*1/stepTDS)/(1/stepTDS);
                    
                    if Iscmax_primer(1)>20*Ip(1) && Iscmax_backup(2)>20*Ip(2)
                        Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
                        Particle(i).TimeOperationSek(2)=Particle(i).TimeOperationPrim(2);
                    else if Iscmax_primer(1)>20*Ip(1) && Iscmax_backup(2)<20*Ip(2)
                            Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
                            Particle(i).TimeOperationSek(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_backup(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                        else if Iscmax_primer(1)<20*Ip(1) && Iscmax_backup(2)<20*Ip(2)
                                Particle(i).TimeOperationPrim(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_primer(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                                Particle(i).TimeOperationSek(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_backup(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                            end
                        end
                    end
                    if Iscmax_primer(2)>20*Ip(2)
                        Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
                    else if Iscmax_primer(2)<20*Ip(2)
                            Particle(i).TimeOperationPrim(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_primer(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                        end
                    end
                else break
                end
            end
            
            while Particle(i).TimeOperationPrim(2)>Target_TopMax
                if TDSmax>=TDSmin
                    TDSmax=TDSmax-stepTDS;
                    Particle(i).xTDS(2)=ceil((((TDSmax-TDSmin)*rand)+TDSmin)*1/stepTDS)/(1/stepTDS);
                    
                    if Iscmax_primer(1)>20*Ip(1) && Iscmax_backup(2)>20*Ip(2)
                        Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
                        Particle(i).TimeOperationSek(2)=Particle(i).TimeOperationPrim(2);
                    else if Iscmax_primer(1)>20*Ip(1) && Iscmax_backup(2)<20*Ip(2)
                            Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
                            Particle(i).TimeOperationSek(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_backup(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                        else if Iscmax_primer(1)<20*Ip(1) && Iscmax_backup(2)<20*Ip(2)
                                Particle(i).TimeOperationPrim(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_primer(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                                Particle(i).TimeOperationSek(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_backup(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                            end
                        end
                    end
                    if Iscmax_primer(2)>20*Ip(2)
                        Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
                    else if Iscmax_primer(2)<20*Ip(2)
                            Particle(i).TimeOperationPrim(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_primer(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                        end
                    end
                else break
                end
            end
        end
        
        for i=1:nPop
            Particle(i).vTDS=0.1*Particle(i).xTDS;
        end
        
        % Error CTI
        for i=1:nPop
            CTI(i)=Particle(i).TimeOperationSek(relay_primary_backup(1,2))-Particle(i).TimeOperationPrim(relay_primary_backup(1,1));
        end
        
        % Objective Function
        for i=1:nPop
            CostFunction(i)=sum(Particle(i).TimeOperationPrim);
            Particle(i).Cost=CostFunction(i);
        end
        
        % Update The Personal Best
        for i=1:nPop
            for j=1:nRelay
                Particle(i).Best.xTDS(j)=Particle(i).xTDS(j);
                Particle(i).Best.TimeOperationPrim(j)=Particle(i).TimeOperationPrim(j);
                Particle(i).Best.TimeOperationSek(j)=Particle(i).TimeOperationSek(j);
            end
            
            %Initial Particle best Cost
            Particle(i).Best.Cost=100;
            % Particle(i).Best.Cost=Particle(i).Cost;
        end
        
        % Update The Global Best
        for i=1:nPop
            if Particle(i).Best.Cost<GlobalBest.Cost
                GlobalBest=Particle(i).Best;
            end
        end
        
        %Lineary Increased c3
        for it = 2 : MaxIt
            coeficient3(it).coef3=coeficient3(it-1).coef3+selisih;
        end
        
        
        % Array to Hold Initial Random Solution of Velocity
        vTDS_initial=[];
        for i=1:nPop
            for j=1:nRelay
                vTDS_initial(i,j)=Particle(i).vTDS(j);
            end
        end
        
        % Array to Hold Best Cost Value on Each Iteration
        BestCosts=zeros(MaxIt,1);
        
        %% =============== Main Loop of PSO ===============%%
        TDSmax=12.5;
        TDSmin=0.1;
        
        for it=1:MaxIt
            for i=1:nPop
                for j=2:nRelay
                    % Update c3
                    c3(it)=coeficient3(it).coef3;
                    
                    % Update R
                    R=Particle(randi([1 nPop],1,1)).xTDS;
                    
                    % Update Damping Coefficients
                    w=wmax-((wmax-wmin).*(it/MaxIt));
                    
                    % Update vTDS
                    Particle(i).vTDS(j)=w*Particle(i).vTDS(j)...
                        +c1*rand.*(Particle(i).Best.xTDS(j) - Particle(i).xTDS(j))...
                        +c2*rand.*(GlobalBest.xTDS(j) - Particle(i).xTDS(j))...
                        +c3(it).*rand.*(R(j)-Particle(i).xTDS(j));
                    
                    % Apply vTDS Limits
                    Particle(i).vTDS=max(Particle(i).vTDS,MinvTDS);
                    Particle(i).vTDS=min(Particle(i).vTDS,MaxvTDS);
                    
                    % Update xTDS with STEP
                    Particle(i).xTDS(j)=ceil((Particle(i).xTDS(j)+Particle(i).vTDS(j))*1/stepTDS)/(1/stepTDS);
                    xTDS(i,j)=Particle(i).xTDS(j);
                    
                    % Apply Lower and Upper Bound TDS Limits
                    Particle(i).xTDS(j)=max(Particle(i).xTDS(j),TDSmin);
                    Particle(i).xTDS(j)=min(Particle(i).xTDS(j),TDSmax);
                    
                    % Evaluation TOP
                    % Evaluasi Top Relay Pertama
                    if Iscmax_primer(1)>20*Ip(1) && Iscmax_backup(2)>20*Ip(2)
                        Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
                        Particle(i).TimeOperationSek(2)=Particle(i).TimeOperationPrim(2);
                    else if Iscmax_primer(1)>20*Ip(1) && Iscmax_backup(2)<20*Ip(2)
                            Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
                            Particle(i).TimeOperationSek(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_backup(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                        else if Iscmax_primer(1)<20*Ip(1) && Iscmax_backup(2)<20*Ip(2)
                                Particle(i).TimeOperationPrim(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_primer(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                                Particle(i).TimeOperationSek(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_backup(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                            end
                        end
                    end
                    % Evaluasi Top Relay Kedua
                    if Iscmax_primer(2)>20*Ip(2)
                        Particle(i).TimeOperationPrim(2)=ConstantSat(2)*Particle(i).xTDS(2);
                    else if Iscmax_primer(2)<20*Ip(2)
                            Particle(i).TimeOperationPrim(2)=((coef_k(2)/coef_beta(2))*(Particle(i).xTDS(2)))/(((Iscmax_primer(1,2)/(Ip(1,2)))^coef_alpha(2))-1);
                        end
                    end
                    
                    % Evaluation Ofun
                    CostFunction(i)=sum(Particle(i).TimeOperationPrim);
                    Particle(i).Cost=CostFunction(i);
                end
            end
            
            % Error CTI
            for i=1:nPop
                CTI(i)=Particle(i).TimeOperationSek(relay_primary_backup(1,2))-Particle(i).TimeOperationPrim(relay_primary_backup(1,1));
            end
            
            % Error Cost Evaluation
            for i=1:nPop
                if Particle(i).TimeOperationPrim(2)<Target_TopMin || Particle(i).TimeOperationPrim(2)>Target_TopMax
                    Particle(i).Cost=100;
                end
            end
            
            for i=1:nPop
                if CTI(i)<Target_CTI
                    Particle(i).Cost=100;
                end
            end
            
            % Update Personal & Global Best
            for i=1:nPop
                
                % Update Personal Best
                if Particle(i).Cost<Particle(i).Best.Cost
                    Particle(i).Best.xTDS=Particle(i).xTDS;
                    Particle(i).Best.TimeOperationPrim=Particle(i).TimeOperationPrim;
                    Particle(i).Best.TimeOperationSek=Particle(i).TimeOperationSek;
                    Particle(i).Best.Cost=sum(Particle(i).Best.TimeOperationPrim);
                    
                end
            end
            
            for i=1:nPop
                % Update Global Best
                if Particle(i).Best.Cost<GlobalBest.Cost
                    GlobalBest=Particle(i).Best;
                    GlobalBest.TimeOperationSek(1)=0;
                end
            end
            
            % Store The Best Cost Value
            BestCosts(it)=GlobalBest.Cost;
            
            % Saving The Relays Data
            for j=1:nRelay
                SaveTopPrim(it,part)=GlobalBest.TimeOperationPrim(j);
                SaveTDS(it,part)=GlobalBest.xTDS(j);
                SaveIscmaxPrim(1,part)=Iscmax_prim(2);
                SaveIscmaxPrimer(1,part)=Iscmax_primer(2);
                SaveIscmaxPrimPu(1,part)=Iscmax_primer(2)/SaveIscmaxPrim(1);
                SavekV(1,part)=kV(2);
                SaveFLA(1,part)=FLA(2);
                SaveRelayPrim(1,part)=part;
                SaveRelayBack(1,part-1)=part;
                SaveRelayBack(1,maxpart)=0;
                SaveCT(1,part)=CT(2);
                SaveIscmaxBack(1,part)=Iscmax_back(2);
                SaveIscmaxBackPu(1,part)=Iscmax_backup(2)/SaveIscmaxPrim(1);
                SaveCurveType(1,part)=CurveType(2);
                SaveCoef_k(1,part)=coef_k(2);
                SaveCoef_alpha(1,part)=coef_alpha(2);
                SaveCoef_beta(1,part)=coef_beta(2);
                SaveConstantSat(1,part)=ConstantSat(2);
            end
            
            for j=1:main_backup_pair
                SaveTopSek(it,part-1)=GlobalBest.TimeOperationSek(j+1);
                SaveTopSek(1,maxpart)=0;
            end
            
            % Saving The Objective Function
            TopConverg(it,1)=sum(SaveTopPrim(it,:));
            
            TDSmax=12.5;
            TDSmin=0.1;
            
            for i=1:nPop
                Particle_TDS_Movement((i+((part-1)*nPop)),it)=Particle(i).xTDS(2);
                Particle_Velocity_Movement((i+((part-1)*nPop)),it)=Particle(i).vTDS(2);
                Particle_TopPrim_Movement((i+((part-1)*nPop)),it)=Particle(i).TimeOperationPrim(2);
                Particle_TopBack_Movement((i+((part-1)*nPop)),it)=Particle(i).TimeOperationSek(2);
            end
            
        end
        
        % warningstatus=zeros(1,1);
        if GlobalBest.TimeOperationPrim(2)>=0.9 || warningstatus==1
            fprintf('=================================\n');
            fprintf('            WARNING!!!\n');
            fprintf('=================================\n');
            fprintf('Waktu Operasi Rele > 0,9 detik\n');
            fprintf('Program Akan Mencari Jenis Kurva Lain\n');
            fprintf('Setuju = 1, Tidak = 0\n');
            warningstatus=input('Masukkan Pilihan : ');
            
            if warningstatus==0
                break
            end
            BestCurveTop(1,CurveType)=GlobalBest.TimeOperationPrim(2);
            BestCurveTDS(1,CurveType)=GlobalBest.xTDS(2);
            if CurveType(2)==1
                CurveType(2)=CurveType(2)+1;
            else if CurveType(2)==2
                    CurveType(2)=CurveType(2)+1;
                else if CurveType(2)==3
                        CurveType(2)=CurveType(2)+1;
                    else if CurveType(2)==4
                            CurveType(2)=CurveType(2)+1;
                        else if CurveType(2)==5
                                CurveType(2)=1;
                            end
                        end
                    end
                end
            end
            
            if CurveType(2)==1
                coef_k(2)=0.14;
                coef_alpha(2)=0.02;
                coef_beta(2)=2.97;
            else if CurveType(2)==2
                    coef_k(2)=13.5;
                    coef_alpha(2)=1;
                    coef_beta(2)=1.5;
                else if CurveType(2)==3
                        coef_k(2)=120;
                        coef_alpha(2)=1;
                        coef_beta(2)=13.33;
                    else if CurveType(2)==4
                            coef_k(2)=80;
                            coef_alpha(2)=2;
                            coef_beta(2)=0.808;
                        else if CurveType(2)==5
                                coef_k(2)=315.2;
                                coef_alpha(2)=2.5;
                                coef_beta(2)=1;
                            end
                        end
                    end
                end
            end
            
            if BestCurveTop~=0
                BestCurveTopMin=min(BestCurveTop);
                for zz=1:5
                    if BestCurveTop(zz)<=BestCurveTopMin
                        xx=zz;
                        GlobalBest.TimeOperationPrim(2)=BestCurveTop(xx);
                        GlobalBest.xTDS(2)=BestCurveTDS(xx);
                        CurveType(2)=xx;
                        
                    end
                end
                
                fprintf('Jenis Kurva Akan Diubah Menjadi Jenis Nomor (%.0f) \n',CurveType(2));
                warningstatus=1;
                
                if CurveType(2)==1
                    coef_k(2)=0.14;
                    coef_alpha(2)=0.02;
                    coef_beta(2)=2.97;
                    ConstantSat(2)=((coef_k(2)/coef_beta(2))*1)/((((ipsat/1))^coef_alpha(2))-1);
                else if CurveType(2)==2
                        coef_k(2)=13.5;
                        coef_alpha(2)=1;
                        coef_beta(2)=1.5;
                        ConstantSat(2)=((coef_k(2)/coef_beta(2))*1)/((((ipsat/1))^coef_alpha(2))-1);
                    else if CurveType(2)==3
                            coef_k(2)=120;
                            coef_alpha(2)=1;
                            coef_beta(2)=13.33;
                            ConstantSat(2)=((coef_k(2)/coef_beta(2))*1)/((((ipsat/1))^coef_alpha(2))-1);
                        else if CurveType(2)==4
                                coef_k(2)=80;
                                coef_alpha(2)=2;
                                coef_beta(2)=0.808;
                                ConstantSat(2)=((coef_k(2)/coef_beta(2))*1)/((((ipsat/1))^coef_alpha(2))-1);
                            else if CurveType(2)==5
                                    coef_k(2)=315.2;
                                    coef_alpha(2)=2.5;
                                    coef_beta(2)=1;
                                    ConstantSat(2)=((coef_k(2)/coef_beta(2))*1)/((((ipsat/1))^coef_alpha(2))-1);
                                end
                            end
                        end
                    end
                end
                warningstatus=0;
                
                SaveTopPrim(MaxIt,part)=GlobalBest.TimeOperationPrim(2);
                SaveTDS(MaxIt,part)=GlobalBest.xTDS(2);
                SaveIscmaxPrim(1,part)=Iscmax_prim(2);
                SaveIscmaxPrimer(1,part)=Iscmax_primer(2);
                SaveIscmaxPrimPu(1,part)=Iscmax_primer(2)/SaveIscmaxPrim(1);
                SavekV(1,part)=kV(2);
                SaveFLA(1,part)=FLA(2);
                SaveRelayPrim(1,part)=part;
                SaveRelayBack(1,part-1)=part;
                SaveRelayBack(1,maxpart)=0;
                SaveCT(1,part)=CT(2);
                SaveIscmaxBack(1,part)=Iscmax_back(2);
                SaveIscmaxBackPu(1,part)=Iscmax_backup(2)/SaveIscmaxPrim(1);
                SaveCurveType(1,part)=CurveType(2);
                SaveCoef_k(1,part)=coef_k(2);
                SaveCoef_alpha(1,part)=coef_alpha(2);
                SaveCoef_beta(1,part)=coef_beta(2);
                SaveConstantSat(1,part)=ConstantSat(2);
                break
            else if warningstatus==0
                    break
                end
            end
        else break
        end
    end
    
    % Saving The Relays Data
    for j=1:main_backup_pair
        dt(j)=GlobalBest.TimeOperationSek(relay_primary_backup(j,2))-GlobalBest.TimeOperationPrim(relay_primary_backup(j,1));
        SaveTargetCTI(1,part-1)=Target_CTI;
        SaveTargetCTI(1,maxpart)=0;
        SaveCTI(1,part-1)=dt(j);
        SaveCTI(1,maxpart)=0;
        SaveErrorCTI(1,part-1)=((SaveCTI(part-1)-SaveTargetCTI(part-1))/SaveTargetCTI(part-1))*100;
        SaveErrorCTI(1,maxpart)=0;
        
    end
    
    SaveIscmaxBack2(1,part-1)=SaveIscmaxBack(part);
    SaveIscmaxBack2(1,maxpart)=0;

    %% ===============Plot TCC=============== %%
    % Plot Kurva TCC Saturasi
    figure;
    iscx1=[(Ip(1)):1:(20*Ip(1))];
    iscx2=[(Ip(2)):1:(20*Ip(2))];
    y1=(coef_k(1)*((GlobalBest.xTDS(1)/coef_beta(1))))./(((iscx1./(Ip(1))).^coef_alpha(1))-1);
    y2=(coef_k(2)*((GlobalBest.xTDS(2)/coef_beta(2))))./(((iscx2./(Ip(2))).^coef_alpha(2))-1);
    kurvatcc=loglog(iscx1,y1,iscx2,y2);
    set(kurvatcc,'linewidth',1);
    hold on;
    kurvaiscprim1=loglog([Iscmax_primer(1) Iscmax_primer(1)],[.01 1e3],'--m');
    set(kurvaiscprim1,'linewidth',1.5);
    hold on;
    kurvaiscsek1=loglog([Iscmax_backup(2) Iscmax_backup(2)],[.01 1e3],'--g');
    set(kurvaiscsek1,'linewidth',1.5);
    hold on;
    kurvasat=loglog([20*Ip(1) (max(SaveIscmaxPrimer)+20000)],[((GlobalBest.xTDS(1))*ConstantSat(1)) ((GlobalBest.xTDS(1))*ConstantSat(1))],[20*Ip(2) (max(SaveIscmaxPrimer)+20000)],[((GlobalBest.xTDS(2))*ConstantSat(2)) ((GlobalBest.xTDS(2))*ConstantSat(2))]);
    set(kurvasat,'linewidth',1);
    hold on
    grid on
    
    legend('Kurva Relay Primer','Kurva Relay Backup','Isc Primer 1','Isc Backup 1');
    set(kurvatcc(1),'Color','red');
    set(kurvasat(1),'Color','red');
    set(kurvatcc(2),'Color','blue');
    set(kurvasat(2),'Color','blue');
    
    xlim([0,max((SaveIscmaxPrimer)+20000)]);
    ylim([0,1e3]);
    title('Kurva TCC');
    xlabel('Arus (A) @base kV');
    ylabel('Waktu (s)');
    
    ax=gca;
    ax.XTick=[0.01 300 500 1000 3000 5000 10000 30000 50000 100000];
    ax.YTick=[0.1 0.3 0.5 1 3 5 10 30 50 100 300 500 1000];
    
    
    %% ===============Saving The TDS=============== %%
    for i=1:nPop
        Particle(i).xTDS(1)=GlobalBest.xTDS(2);
        Particle(i).TimeOperationPrim(1)=GlobalBest.TimeOperationPrim(2);
    end
    
    for j=1:main_backup_pair
        Target_TopMin=GlobalBest.TimeOperationPrim(j+1);
        Iscmax_prim(1,j)=Iscmax_prim(j+1);
        Iscmax_back(1,j)=Iscmax_back(j+1);
        kV(1,j)=kV(j+1);
        FLA(1,j)=FLA(j+1);
        %         VR(1,j)=VR(j+1);
        CT(1,j)=CT(j+1);
        ConstantSat(j)=ConstantSat(j+1);
        CurveType(j)=CurveType(j+1);
        coef_k(j)=coef_k(j+1);
        coef_alpha(j)=coef_alpha(j+1);
        coef_beta(j)=coef_beta(j+1);
    end
    
    if part==maxpart
        break
    end
    %% ===============Input Relays Data=============== %%
    part=part+1;
    fprintf('\n');
    fprintf('Masukkan Tegangan Relay ke %.0f (kV)',part);
    kV(1,2)=input(' : ');
    fprintf('Masukkan Primer CT ke %.0f ',part);
    CT(1,2)=input(' : ');
    fprintf('Masukkan FLA relay ke %.0f (A)',part);
    FLA(1,2)=input(' : ');
    fprintf('Masukkan Iscmax Relay Primer ke %.0f (A)',part);
    Iscmax_prim(1,2)=input(' : ');
    fprintf('Masukkkan Iscmax Relay Backup ke %.0f (A)',part);
    Iscmax_back(1,2)=input(' : ');
    fprintf('Jenis Kurva :\n');
    fprintf('(1) Standard Inverse\n');
    fprintf('(2) Very Inverse\n');
    fprintf('(3) Long Time Inverse\n');
    fprintf('(4) Extremely Inverse\n');
    fprintf('(5) Ultra Inverse\n');
    CurveType(2)=input('Pilih Jenis Kurva : ');
    
    if CurveType(2)==1
        coef_k(2)=0.14;
        coef_alpha(2)=0.02;
        coef_beta(2)=2.97;
    else if CurveType(2)==2
            coef_k(2)=13.5;
            coef_alpha(2)=1;
            coef_beta(2)=1.5;
        else if CurveType(2)==3
                coef_k(2)=120;
                coef_alpha(2)=1;
                coef_beta(2)=13.33;
            else if CurveType(2)==4
                    coef_k(2)=80;
                    coef_alpha(2)=2;
                    coef_beta(2)=0.808;
                else if CurveType(2)==5
                        coef_k(2)=315.2;
                        coef_alpha(2)=2.5;
                        coef_beta(2)=1;
                    end
                end
            end
        end
    end

    if kV(2)~=kV(1)
        fprintf('==========================================\n');
        fprintf('                PERHATIAN                  \n');
        fprintf('==========================================\n');
        fprintf('LEVEL TEGANGAN BERBEDA TELAH TERDETEKSI \n');
        fprintf('Apakah CTI Ingin Diminimalkan?\n');
        fprintf('YA = 1, TIDAK = 2\n');
        CTIchoice=input('Masukkan Pilihan : ');
        if CTIchoice==1
            fprintf('STATUS : CTI DIMINIMALKAN\n');
            fprintf('==========================================');
            Target_TopMin=0.1;
            Target_CTI=0;
            %             Target_CTI=0.01;
            fprintf('\n');
        else if CTIchoice==2
                Target_CTI=input('Target CTI (s) : ');
                Target_TopMin=0.1;
                fprintf('\n');
            end
        end
    else Target_CTI=input('Target CTI (s) : ');
    end
    
    relay_primary_backup=[1 2];
    
    TDSmax=12.5;
    TDSmin=0.1;
    

end
%% ===============Results=============== %%
fprintf('\n');
fprintf('\n');
disp('=======================================================================');
disp('                               Results                           ');
disp('=======================================================================');
for it=1:MaxIt
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(TopConverg(it))]);
end
disp('Sum Top Primer :');
disp(TopConverg(MaxIt,:));
fprintf('\n');
fprintf('\n')


disp('==========================================================================================');
disp('                                 Relays Data Summary                           ');
disp('==========================================================================================');
disp('Relay     Voltage(kV)    Primary CT      Ip (A)      Tap        TDS         Curve Type');
for j=1:maxpart
    fprintf('%3.0f',SaveRelayPrim(j));
    fprintf('%15.2f',SavekV(j));
    fprintf('%13.0f',SaveCT(j));
    fprintf('%16.2f',1.05*SaveFLA(j));
    fprintf('%10.2f',(1.05*SaveFLA(j)/SaveCT(j)));
    fprintf('%11.2f',SaveTDS(MaxIt,j));
    if SaveCurveType(j)==1
        fprintf('      Standard Inverse');
    else if SaveCurveType(j)==2
            fprintf('      Very Inverse');
        else if SaveCurveType(j)==3
                fprintf('      Long Time Inverse');
            else if SaveCurveType(j)==4
                    fprintf('      Extremely Inverse');
                else if SaveCurveType(j)==5
                        fprintf('      Ultra Inverse');
                    end
                end
            end
        end
    end
    
    fprintf('\n');
end
fprintf('\n');
fprintf('\n');

disp('==========================================================================================================');
disp('                                         Time Operation Analysis                            ');
disp('==========================================================================================================');
disp('Relay Primer    Relay Backup   Iscmax Primer(A)  Iscmax Backup(A)   Top Primer(s)   Top Backup(s)   CTI(s)');
for j=1:maxpart
    fprintf('%6.0f',SaveRelayPrim(j));
    fprintf('%17.0f',SaveRelayBack(j));
    fprintf('%18.2f',SaveIscmaxPrim(j));
    fprintf('%18.2f',SaveIscmaxBack2(j));
    fprintf('%18.4f',SaveTopPrim(MaxIt,j));
    fprintf('%15.4f',SaveTopSek(MaxIt,j));
    fprintf('%14.4f',SaveCTI(j));
    fprintf('\n');
end

fprintf('\n');
fprintf('\n');

disp('=================================================================');
disp('                            CTI Analysis                   ');
disp('=================================================================');
disp('Relay Primer   Relay Backup   CTI(s)   Target CTI   Difference(%)');
for j=1:maxpart-1
    fprintf('%6.0f',SaveRelayPrim(j));
    fprintf('%15.0f',SaveRelayBack(j));
    fprintf('%15.4f',SaveCTI(j));
    fprintf('%12.4f',SaveTargetCTI(j));
    fprintf('%11.4f',SaveErrorCTI(j));
    fprintf('\n');
end


figure;
plot(TopConverg, 'LineWidth',2);
xlabel('Iterations');
ylabel('Best Costs');
grid on
title('Kurva Konvergensi');

figure;
for j=1:maxpart
    iscx2=[(1.05*SaveFLA(j)*SavekV(j)/kVbase):1:(20*SaveFLA(j)*SavekV(j)/kVbase)];
    y2=(SaveCoef_k(j)*((SaveTDS(MaxIt,j)/SaveCoef_beta(j))))./(((iscx2./(1.05*SaveFLA(j)*SavekV(j)/kVbase)).^SaveCoef_alpha(j))-1);
    kurvatcctot2=loglog(iscx2,y2);
    hold on
    grid on
end

xlim([0,max((SaveIscmaxPrimer)+20000)]);
ylim([0,1e3]);
title('Kurva TCC');
xlabel('Arus (A) @base kV');
ylabel('Waktu (s)');
ax=gca;
ax.XTick=[0.01 300 500 1000 3000 5000 10000 30000 50000 100000];
ax.YTick=[0.1 0.3 0.5 1 3 5 10 30 50 100 300 500 1000];

for j=1:maxpart
    PlotLegendRelay{j}=strcat('Kurva Relay  ', num2str(j));
end

title('Kurva TCC Logaritmik');
xlabel('Arus (A)');
ylabel('Waktu (s)');

