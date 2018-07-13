DataFilePath = 'C:/Users/KevinYeung/Desktop/intern Mont/examples_cvs/';
addpath('C:/Users/KevinYeung/Desktop/intern Mont/examples_cvs/');
data_list=struct2cell(dir(fullfile(DataFilePath)));
file_name_list = data_list(1,3:end);
% DataFileName = '20150210_100X_Movie_Calibration_RI510_128_noT_001_visit_2_c2_CALIB_RF.csv';
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
% DureeSimu = 3000; % 1200 secondes Attention cette valeur ne doit pas ?tre inf?rieure a la dur?e exp?rimentale
%CoeffNormalisation = 1
Gap4TrainDef = 2;

exp_data = {};
simu_data = {};

for data_i = 1:8%length(file_name_list) 
    % -------input file----------
    DataFileName = file_name_list{data_i}
    DataExp = table2array(readtable(DataFileName, 'Delimiter','semi'));  
    DataExpSmooth = transp(DataExp(:,2));  % signal intensity in experiment data
    exp_data{end+1} = DataExpSmooth;
    % -------Parameters----------
    DureeSimu = DataExp(end,1)+DataExp(2,1); % (s)
    frame_num = length(DataExp(:,1));
    DureeSignal = round((TaillePreMarq+TailleSeqMarq + TaillePostMarq) / Polym_speed); % (s)
    DureeAnalysee = DureeSignal + DureeSimu - round(TaillePreMarq / Polym_speed); % (s)
    FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed); % how many interval(possible poly start position) in 1s
    num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed));
    Parameters = {FreqEchSimu, FreqEchImg,DureeAnalysee,NbrSondeFluo,...
            ProbeByIntensitie_nb,TaillePreMarq,TailleSeqMarq, TaillePostMarq, Polym_speed,frame_num};
    % ---------estimate poly number---------
    eval_poly_num = [];
    GA_fitness_try_poly = @(x) sum((sumSignal(x,Parameters)-DataExpSmooth).^2);
    for num_p = 30:30:900
        fitmin = 1e+10;
        Estimate_polys = zeros(50,num_p);
        for i = 1:50
            poly_posi = randperm(num_possible_poly,num_p); % randomly choose 25 poly position (test1)
            fit_now = GA_fitness_try_poly(poly_posi);
            if fit_now<fitmin
                fitmin = fit_now;
            end
        %     Trans_positions_art = find(Pattern_poly_art(i,:)==1);
        end
        eval_poly_num(end+1) = fitmin;  
    end
    Nbr_poly_estimate = find(eval_poly_num==min(eval_poly_num))*30

    % -----------initial population for GA----------
    Nbr_simu_DNA = 500; % number of "chromosome" in GA
%     Nbr_poly_estimate = 270; % different for every data
    Pattern_polys = zeros(Nbr_simu_DNA,num_possible_poly);
    for i = 1:Nbr_simu_DNA
        Pattern_polys(i,randperm(num_possible_poly,Nbr_poly_estimate)) = 1; % randomly choose 25 poly position (test1)
    %     Trans_positions_art = find(Pattern_poly_art(i,:)==1);
    end
    
    % ----------GA----------------
    GA_fitness_exp = @(x) sum((sumSignal(find(x==1),Parameters)-DataExpSmooth).^2);
    options = gaoptimset;
    options = gaoptimset(options,'PopulationType', 'bitstring', 'CreationFcn',@gacreationuniform,'MutationFcn',@mutationuniform,...
        'Display','iter','TolFun',1e-8,'Paretofraction',0.35,'Generations',200);

    options = gaoptimset(options,'PopulationSize', Nbr_simu_DNA);
    options = gaoptimset(options,'InitialPopulation', Pattern_polys);
    % x_GA = ga(GA_fitness_art,length(Pattern_poly_art),[],[],[],[],[],[],[],[],options);
    x_GA_exp = ga(GA_fitness_exp,length(Pattern_polys(1,:)),[],[],[],[],[],[],[],[],options);
    simu_data{end+1} = x_GA_exp;
    save(['simu_data(',num2str(data_i),')_200iters_500population'],'x_GA_exp'); 
end

%% visualize
for ii = 5:8
   figure(ii)
   plot(sumSignal(find(simu_data{ii}==1),Parameters))
   hold on 
   plot(exp_data{ii},'red')
end

%% Add noise
% study the noise from experiment data and simulation data
% each row of signal_simu is the sumSignal of one simu data of
% corresponding exp data (1 file)
signal_simu = sumSignal(find(simu_data{1}==1),Parameters); 
for ii = 2:length(simu_data) 
    signal_simu(ii,:) = sumSignal(find(simu_data{ii}==1),Parameters);
end 
signal_exp = cell2mat(transpose(exp_data)); % each row is one exp file data
diff_simu_exp = signal_simu-signal_exp;
var_window = 5;
varDiff_simuSignal = zeros(2,length(round(min(signal_simu(:))):var_window:round(max(signal_simu(:))))); % store var of 1 window and signal intensity of this window
idx = 1;
for ii = round(min(signal_simu(:))):var_window:round(max(signal_simu(:)))
    df_in_window = diff_simu_exp(find(signal_simu>ii & signal_simu<ii+var_window));
    varDiff_simuSignal(1,idx) = ii+var_window/2;
    varDiff_simuSignal(2,idx) = var(df_in_window);
    idx = idx+1;
end

% -----visualize noise in terms of signal intensity---------
figure(3)
plot(varDiff_simuSignal(1,:),varDiff_simuSignal(2,:),'o')
plot(signal_simu,(signal_simu-signal_exp),'o')

% -------regression----------
x = varDiff_simuSignal(1,:);
y = varDiff_simuSignal(2,:);
p = polyfit(x,y,1);
f = polyval(p,x); 
figure(4)
plot(x,y,'o',x,f,'x') 

% --------add noise to simu data-------
mean_sig = mean(diff_simu_exp(:));
var_sig = polyval(p,signal_simu(4,:));
simu_signal_noise = signal_simu(4,:) + mean_sig + (2*(rand(1,length(signal_simu(4,:)))-0.5)).*var_sig;
figure(5)
plot(simu_signal_noise)
hold on 
plot(signal_exp(4,:),'red')

sum((simu_signal_noise-exp_data{4}).^2)
sum((signal_simu(4,:)-exp_data{4}).^2)

% -----create artificial data and add noise---------
Nbr_poly = 90;
Pattern_poly_art = zeros(1,num_possible_poly);
Pattern_poly_art(randi(num_possible_poly,1,Nbr_poly)) = 1; % randomly choose 25 poly position (test1)
Trans_positions_art = find(Pattern_poly_art==1);
sum_signal_art = sumSignal(Trans_positions_art,Parameters);
% add noise
var_sig = polyval(p,sum_signal_art);
% amplitude_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1]; % control amplitude of noise
amplitude_list = [0];
% try artificial data for 5 different amplitude value
art_signal_noise = zeros(length(amplitude_list),length(sum_signal_art));
for ii = 1:length(amplitude_list)
    art_signal_noise(ii,:) = sum_signal_art + mean_sig + (amplitude_list(ii)*2*(rand(1,length(sum_signal_art))-0.5)).*var_sig;
end

%% GA again to estimate artificial data
% -----------initial population for GA----------
% try decovolute artificial data for 5 different amplitude value
Pattern_poly_simu_diff_polyNbr = zeros(10,num_possible_poly);
Pattern_poly_art_diff_polyNbr = zeros(10,num_possible_poly);
art_signal_noise_diff_polyNbr = zeros(10,frame_num);
for ii = 1:10
    %create artificial data for different poly nbr
    Nbr_poly = ii*30;
    Pattern_poly_art = zeros(1,num_possible_poly);
    Pattern_poly_art(randperm(num_possible_poly,Nbr_poly)) = 1; % randomly choose 25 poly position (test1)
    Trans_positions_art = find(Pattern_poly_art==1);
    sum_signal_art = sumSignal(Trans_positions_art,Parameters);
    
    Pattern_poly_art_diff_polyNbr(ii,:) = Pattern_poly_art;  % store art_data poly position

    % add noise
    var_sig = polyval(p,sum_signal_art);
    % amplitude_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1]; % control amplitude of noise
    amplitude = 0.5;
    % try artificial data for different amplitude value
    art_signal_noise_i = sum_signal_art + mean_sig + (amplitude*2*(rand(1,length(sum_signal_art))-0.5)).*var_sig;
    art_signal_noise_diff_polyNbr(ii,:) = art_signal_noise_i;
    % prepare initial population 
    Nbr_simu_DNA = 500; % number of "chromosome" in GA
    Nbr_poly_estimate = ii*30-10; % different for every data
    Pattern_polys = zeros(Nbr_simu_DNA,num_possible_poly);
    for i = 1:Nbr_simu_DNA
        Pattern_polys(i,randperm(num_possible_poly,Nbr_poly_estimate)) = 1; % randomly choose poly position 
    %     Trans_positions_art = find(Pattern_poly_art(i,:)==1);
    end

    % ----------GA----------------
    GA_fitness_art_noise = @(x) sum((sumSignal(find(x==1),Parameters)-art_signal_noise_i).^2);
%     GA_fitness_art = @(x) sum((sumSignal(find(x==1),Parameters)-sum_signal_art).^2);
    options = gaoptimset;
    options = gaoptimset(options,'PopulationType', 'bitstring', 'CreationFcn',@gacreationuniform,'MutationFcn',@mutationuniform,...
        'Display','iter','TolFun',1e-8,'Paretofraction',0.35,'Generations',300);

    options = gaoptimset(options,'PopulationSize', Nbr_simu_DNA);
    options = gaoptimset(options,'InitialPopulation', Pattern_polys);
    % x_GA = ga(GA_fitness_art,length(Pattern_poly_art),[],[],[],[],[],[],[],[],options);
    x_GA_art_noise = ga(GA_fitness_art_noise,length(Pattern_polys(1,:)),[],[],[],[],[],[],[],[],options);
    Pattern_poly_simu_diff_polyNbr(ii,:) = x_GA_art_noise;
    save(['simu_art_data_noise_diff_polyNbr(',num2str(ii*30),')_300iters_500population'],'x_GA_art_noise'); 
end
% -----visualize GA result------ 
% blue GA data, red artifical experiment data
figure(1)
for ii = 1:10
    subplot(4,3,ii)
    plot(sumSignal(find(Pattern_poly_art_diff_polyNbr(ii,:)==1),Parameters),'red')
    hold on 
    plot(sumSignal(find(Pattern_poly_simu_diff_polyNbr(ii,:)==1),Parameters))
    plot(art_signal_noise_diff_polyNbr(ii,:),'black')
    legend('art data','simu data','art data noise')
    title(['Nbr Poly=',num2str(ii*30),', amp=0.5'])
end 


Pattern_poly_art_no_noise = x_GA_art;
Pattern_poly_list_diff_amp = Pattern_poly_list; % store data for different amp value, load data from file
Pattern_poly_art_diff_amp = Pattern_poly_art; % Pattern poly for num_poly = 90, load data from file
distance_pattern_poly = @(x) sum((x(1:length(Trans_positions_art))-Trans_positions_art).^2);
figure(1)
plot(find(Pattern_poly_list(1,:)==1),1,'Marker','o','color','red')
hold on 
plot(Trans_positions_art,1,'Marker','x','color','blue')

% ----------Measure the difference between art and simu data---------
% use Hausdorff distance

% different poly nbr
hd_dist_list = [];
df_nbr_poly = [];
for ii = 1:10 
    art_posi = find(Pattern_poly_art_diff_polyNbr(ii,:));
    simu_posi = find(Pattern_poly_simu_diff_polyNbr(ii,:));
    hd_dist_list(ii) = HausdorffDist(transpose(art_posi),transpose(simu_posi));
    df_nbr_poly(ii) = abs(length(art_posi)-length(simu_posi))/length(art_posi);   % difference of poly number
end

plot(hd_dist_list,'o')
hold on 
plot(df_nbr_poly*200,'marker','o','color','red')

% different amplitude









    
