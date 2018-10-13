%% Import files and read experiment data
%DataFilePath = 'C:/Users/KevinYeung/Desktop/intern Mont/examples_cvs/';
DataFileName = '20150210_100X_Movie_Calibration_RI510_128_noT_001_visit_2_c2_CALIB_RF.csv';
%cd(DataFilePath);
DataExp = table2array(readtable(DataFileName, 'Delimiter','semi'));  
DataExpSmooth = transp(DataExp(:,2));  % signal intensity in experiment data

%ScripPath = 'C:\Users\KevinYeung\Desktop\intern Mont';
%cd(ScripPath);
%% Model and Parameters 
%--------transcription--------
NbrSondeFluo = 128; % 128 sondes per mRNA (ref publi)
TaillePreMarq = 700; % 700 bases (ref publi)
TailleSeqMarq = 5800; % 2900 bases (ref publi)
TaillePostMarq = 1600; % 1600 bases (ref publi)
DureeProcessing = 100; %  en secondes (100s temps moyen ref publi)
% --------polymerases------- 
EspaceInterPolyMin = 30; % en base (40 bases)
Polym_speed = 67; % average speed bases par seconde (Ref publi)
ProbeByIntensitie_nb = NbrSondeFluo; % how many probes are needed to observe 1 intensity
% -------sample freqence------
NombreRandomIter = 1000; %at least 1000 but 10000 is fine for 1200 secondes et FreqEchSimu1
FreqEchImg = (1/3); % 1/3 image par seconde
% FreqEchSimu = 2; % nb de calcul par seconde %% Attention cette donn?e doit ?tre enti?re. -> int ; no float ! ET PROP a la dur?e des calculs
DureeSimu = 3000; % 1200 secondes Attention cette valeur ne doit pas ?tre inf?rieure a la dur?e exp?rimentale
%CoeffNormalisation = 1
Gap4TrainDef = 2;
DureeSimu = DataExp(end,1)+DataExp(2,1); % (s)
frame_num = length(DataExp(:,1));
DureeSignal = round((TailleSeqMarq + TaillePostMarq) / Polym_speed); % (s)
DureeAnalysee = DureeSignal + DureeSimu - round(TaillePreMarq / Polym_speed); % (s)
FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed);
Parameters = {FreqEchSimu, FreqEchImg,DureeSimu,NbrSondeFluo,...
            ProbeByIntensitie_nb,TaillePreMarq,TailleSeqMarq, TaillePostMarq, Polym_speed,frame_num};
num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed));

%% estimate how many polymerase on DNA by "counting" peaks of exp. signal
Nbr_poly = 25; % assumption 
% -----create artificial data---------
Pattern_poly_art = zeros(1,num_possible_poly);
Pattern_poly_art(randi(num_possible_poly,1,Nbr_poly)) = 1; % randomly choose 25 poly position (test1)
Trans_positions_art = find(Pattern_poly_art==1);
sum_signal_art = sumSignal(Trans_positions_art,Parameters);

%% Initiation 
% set Nstart DNA, for each DNA there are Nbr_poly polymerase which are
% ready to transcript. Choose the best DNA 
Nbr_simu_DNA = 500; % number of "chromosome" in GA
Nbr_poly_estimate = 25; % test1
Pattern_polys = zeros(Nbr_simu_DNA,num_possible_poly);
for i = 1:Nbr_simu_DNA
    Pattern_polys(i,randi(num_possible_poly,1,Nbr_poly_estimate)) = 1; % randomly choose 25 poly position (test1)
%     Trans_positions_art = find(Pattern_poly_art(i,:)==1);
end
% try a objective function
% objFunGA(Pattern_polys(10,:),sum_signal_art,Parameters)

%% Solve the problem by genetic algo
% start with an easy case
% f = @(x) sumSignal(x,Parameters)-
f = @(x) sum((x-[1,0,0,1]).^2);
% options = gaoptimset;
% options = gaoptimset(options,'PopulationType', 'bitstring');
% options = gaoptimset(options,'PopulationSize', 20);
% options = gaoptimset(options,'InitialPopulation', [1,0,0,1]);
% xt = ga(f,4,[],[],[],[],[],[],[],[],options);
x = lsqnonlin(GD_fitness,[0.5,])

% ------GA------
% GA_fitness = @(x) sum((sumSignal(find(x==1),Parameters)-DataExpSmooth).^2);
GA_fitness_art = @(x) sum((sumSignal(find(x==1),Parameters)-sum_signal_art).^2);
options = gaoptimset;
options = gaoptimset(options,'PopulationType', 'bitstring', 'CreationFcn',@gacreationuniform,'MutationFcn',@mutationuniform,...
    'Display','iter','TolFun',1e-6,'Paretofraction',0.35,'Generations',100);

options = gaoptimset(options,'PopulationSize', Nbr_simu_DNA);
options = gaoptimset(options,'InitialPopulation', Pattern_polys);
x_GA = ga(GA_fitness_art,length(Pattern_poly_art),[],[],[],[],[],[],[],[],options);

% ----Gradient descent-------
GD_fitness = @(x) sumSignal(x,Parameters)-sum_signal_art; % x: positions of poly
trans_position_GD = find(x_GA==1);
op_gd = optimset(@lsqnonlin);
op_gd = optimset(op_gd,'Display','iter','TolFun',1e-1000,'TolX',1e-1000,'MaxIter',500,'MaxFunEvals',8000);
shift_window = 100;
% x_GD = lsqnonlin(GD_fitness,trans_position_GD,trans_position_GD-shift_window,trans_position_GD+shift_window,op_gd);
x_GD = lsqnonlin(GD_fitness,trans_position_GD+10,[],[],op_gd);


% -------------Yang's GD-----------------
% for every transcription start position (x_GA==1), try to shift it (left
% or right), to optimize the result
GD_y_fitness = @(x) sum((sumSignal(x,Parameters)-sum_signal_art).^2); % x: positions of poly
shift_window = 50; % one position can move [-s_w,s_w]
trans_position_y = find(x_GA==1);
Min_fit = GD_y_fitness(trans_position_y);
for posi_i = 1: length(trans_position_y)
    new_tr_p = trans_position_y;
    for j = -shift_window:shift_window
        new_tr_p(posi_i) = trans_position_y(posi_i)+j;
        if new_tr_p(posi_i) <= 0
            continue
        end
        if GD_y_fitness(new_tr_p)<Min_fit
            trans_position_y = new_tr_p;
            Min_fit = GD_y_fitness(trans_position_y);
        end
    end
end
% trans_position =  

% -----visualize GA result------ 
% blue GA data, red artifical experiment data
figure(1)
plot(sumSignal(find(x_GA==1),Parameters))
hold on 
plot(sum_signal_art,'red')
hold on
plot(sumSignal(trans_position_y,Parameters),'green')
figure(2)
plot(find(x_GA==1),1,'Marker','o','color','red')
hold on 
plot(Trans_positions_art,1,'Marker','x','color','blue')

