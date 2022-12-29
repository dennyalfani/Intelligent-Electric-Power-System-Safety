clc
clear all

global iscp iscb ip_base aa bb kk

fprintf('**Perhitungan TDS berbasis Multiobjective GA**\n');
fprintf('pada jenis kurva Extreamly Inverse\n\n');

fprintf('-Nilai kV base yang dipiih 6 kV\n');
fprintf('-Nilai I pick up = 1.05 FLA\n\n');

fprintf('Masukkan data relay\n');
id = input('Masukkan ID Relay : ','s');
v = input('Masukkan Tegangan Operasi : ');
fla = input ('Masukkan FLA : ');
isc_p = input ('Masukkan Isc Primer : ');
isc_b = input ('Masukkan Isc Backup : ');
top = input ('Masukkan target Top Primer : ');
cti = input ('Masukkan nilai CTI : ');
ctp = input ('Masukkan nilai Primer CT : ');

%Parameter Kurva Extreamly Inverse
aa=2; bb=0.808; kk=80;

kvbase=6;

fla_base = fla*v/kvbase;
iscp_base = isc_p*v/kvbase;
iscb_base = isc_b*v/kvbase;
isc_sat = fla_base*20;
ip_base = 1.05*fla_base;

%mencari nilai isc primer terhadap saturasi
if iscp_base > isc_sat
    iscp = isc_sat;
else
    iscp = iscp_base;
end
%iscp

%mencari nilai isc primer terhadap saturasi
if iscb_base > isc_sat
    iscb = isc_sat;
else
    iscb = iscb_base;
end
%iscb


FitnessFunction = @fun;
numberOfVariables = 2;
lb = [(top) 0.1]; % Lower bound
ub = [1 1]; % Upper bound
A = [1,-1]; % No linear inequality constraints
b = [-(cti)]; % No linear inequality constraints
Aeq = []; % No linear equality constraints
beq = []; % No linear equality constraints
options = optimoptions(@gamultiobj,'PlotFcn',@gaplotpareto);

[x,Fval,exitFlag,Output] = gamultiobj(FitnessFunction,numberOfVariables,A, ...
    b,Aeq,beq,lb,ub,options);
disp('Hasil Optimisisasi'),disp(id)
top_optimum_primer = x(:,1)
top_optimum_backup = x(:,2)
time_dial = Fval(:,1)
pick_up = 1.05*fla/ctp

function y = fun (x)

global iscp iscb ip_base aa bb kk

y(1) =(x(1)*((((iscp/ip_base)^aa)-1)*bb))/kk;
y(2) =(x(2)*((((iscb/ip_base)^aa)-1)*bb))/kk;
y(1) = y(2);
y = [y(1) y(2)];
end