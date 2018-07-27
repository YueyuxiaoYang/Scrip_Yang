%% GA again to estimate artificial data
% try decovolute artificial data for different GA 
%iterations[50:50:300], rep=5, amp = 0.5, polyNbr = 150
simu_data_diff_gaIters = struct;
art_data_diff_gaIters = struct;
count_simu_gaIters = 1;
Nbr_poly = 150;
Pattern_poly_art = zeros(1,num_possible_poly);
Pattern_poly_art(randperm(num_possible_poly,Nbr_poly)) = 1; % randomly choose 25 poly position (test1)
Trans_positions_art = find(Pattern_poly_art==1);
sum_signal_art = sumSignal(Trans_positions_art,Parameters);

% add noise
var_sig = polyval(p,sum_signal_art);
% amplitude_list = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1]; % control amplitude of noise
amplitude = 0.5;
% try artificial data for different amplitude value
art_signal_noise_i = sum_signal_art + mean_sig + (amplitude*2*(rand(1,length(sum_signal_art))-0.5)).*var_sig;
%art_signal_noise_diff_polyNbr(ii,:) = art_signal_noise_i;

art_data_diff_gaIters.Amp = 0.5;
art_data_diff_gaIters.polyNbr = 150;
art_data_diff_gaIters.trans_posi_art = Trans_positions_art;
art_data_diff_gaIters.signal_noise = art_signal_noise_i;
save(['art_diff_gaIters(50-300)_pN150_amp05_100iters_500population'],'art_data_diff_gaIters'); 


for ga_i = 1:6
    %create artificial data for different poly nbr
    ga_iters = ga_i*50;
  
    for rep_i = 1:5
        % prepare initial population 
        Nbr_simu_DNA = 500; % number of "chromosome" in GA
        Nbr_poly_estimate = 150-20; % different for every data
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
            'Display','iter','TolFun',1e-10,'Paretofraction',0.35,'Generations',ga_iters);

        options = gaoptimset(options,'PopulationSize', Nbr_simu_DNA);
        options = gaoptimset(options,'InitialPopulation', Pattern_polys);
        % x_GA = ga(GA_fitness_art,length(Pattern_poly_art),[],[],[],[],[],[],[],[],options);
        x_GA_art_noise = ga(GA_fitness_art_noise,length(Pattern_polys(1,:)),[],[],[],[],[],[],[],[],options);
        
        simu_data_diff_gaIters(count_simu_gaIters).polyNbr = Nbr_poly;
        simu_data_diff_gaIters(count_simu_gaIters).amp = 0.5;
        simu_data_diff_gaIters(count_simu_gaIters).trans_posi = find(x_GA_art_noise==1);
        count_simu_gaIters = count_simu_gaIters+1;
        save(['simu_diff_gaIters(',num2str(ga_i),')_pN150_amp05_100iters_500population'],'simu_data_diff_gaIters'); 
%         save(['art_data_noise_diff_polyNbr(polyNbr',num2str(ga_i*30),')_100iters_500population'],'art_data_diff_polyNbr');
    end
end

% ============different GA maxIters==========
jaccard_dist_gaIter = zeros(5,6);
jaccard_dist_gaIter_approx= zeros(5,6);
fit_gaIter = zeros(5,6);
fit_fuc_gaIter = @(x)sum((sumSignal(x,Parameters)-art_signal_noise_i).^2)
jd_dist = @(a,b) 1-length(intersect(a,b))/length(union(a,b)); % compute Jaccard distance
jd_approximate_dist =  @(ref,b) 1-(length(intersect(unique([ref-1,ref-2,ref,ref+1,ref+2]),unique([b-1,b-2,b,b+1,b+2])))/length(union(unique([ref-1,ref-2,ref,ref+1,ref+2]),unique([b-1,b-2,b,b+1,b+2]))));

for ii = 1:length(simu_data_diff_gaIters) 
    s = simu_data_diff_gaIters(ii);
    art_posi = art_data_diff_gaIters.trans_posi_art;
    jaccard_dist_gaIter(ii) = jd_dist(art_posi,s.trans_posi);
    jaccard_dist_gaIter_approx(ii) = jd_approximate_dist(art_posi,s.trans_posi);
    fit_gaIter(ii) = fit_fuc_gaIter(s.trans_posi);
%     s_a = simu_data_diff_amp_rep(ii);
%     art_posi_amp = find(Pattern_poly_art_diff_amp==1);
%     jd_GA_amp(ii) = jd_dist(art_posi_amp,s_a.trans_posi);
%     
end

figure(33)
hold on
boxplot(jaccard_dist_gaIter,'color','r','Labels',50:50:300)
boxplot(jaccard_dist_gaIter_approx,'color','b','Labels',50:50:300)
plot(mean(jaccard_dist_gaIter),'marker','o','color','r')
plot(mean(jaccard_dist_gaIter_approx),'marker','o','color','b')
legend('Jaccard distance after GA','Jaccard approximate distance after GA')
title('Compare Jaccard distance of 6 GA MaxIteration(30-300) with 5 repeat)')
xlabel('Max Iteration of GA') 
ylabel('Jaccard distance')
ylim([0.7,1])

figure(34)
hold on 
boxplot(fit_gaIter,'color','g','Labels',50:50:300)
title('Compare fitness in GA of 6 GA MaxIteration(30-300) with 5 repeat)')
xlabel('Max Iteration of GA') 
ylabel('Fitness Value')