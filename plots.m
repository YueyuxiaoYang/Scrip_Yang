% =================plots for ppt or pdf======================
% visualize the simulation results
% ----------------load data-----------------------------
% 
load('data_simu_727');

%% ----------------- initial events ------------------------
% an example
X = art_data_diff_polyNbr(1).trans_posi_art; % artificial data
Y = 0:0.1:1;
[x,y] = meshgrid(X,Y);
line(x,y,'color','b')
xlabel('time(s)')
ylabel('initial event')

%----------------- signal intensity(sumSignal)--------------
% signal from 1 initial event(getSignal), artificial data
plot(sumSignal(X(3),Parameters))
xlim([0,100])
ylim([0,1.1])
xlabel('time(s)')
ylabel('signal intensity')
% signal of a set of initial events(sumSigal), artificial data
plot(sumSignal(X,Parameters))
xlabel('time(s)')
ylabel('signal intensity')

%% ---------- simulation with different GA maxIter -----------
figure()
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

figure()
hold on 
boxplot(fit_gaIter,'color','g','Labels',50:50:300)
title('Compare fitness in GA of 6 GA MaxIteration(30-300) with 5 repeat)')
xlabel('Max Iteration of GA') 
ylabel('Fitness Value')

%% ----- Jaccard distance for different polymerase number------ 
% use <jaccard_dist_GA_polyNbr_rep> and <jaccard_dist_GA_polyNbr_rep_approx>
figure()
hold on
boxplot(jaccard_dist_GA_polyNbr_rep,'color','r')
boxplot(jaccard_dist_GD_polyNbr_rep,'color','b')
plot(mean(jaccard_dist_GA_polyNbr_rep),'marker','o','color','r')
plot(mean(jaccard_dist_GD_polyNbr_rep),'marker','o','color','b')
legend('Jaccard distance after GA','Jaccard distance after GD')
title('Compare Jaccard distance of 10 PolyNbr(30-300) with 8 repeat)')
xlabel('Polymerase Number(x30)')
ylabel('Jaccard distance')
ylim([0.8,1])

figure()
boxplot(jaccard_dist_GA_polyNbr_rep_approx,'color','r')
hold on
boxplot(jaccard_dist_GD_polyNbr_rep_approx,'color','b')
plot(mean(jaccard_dist_GA_polyNbr_rep_approx),'marker','o','color','r')
plot(mean(jaccard_dist_GD_polyNbr_rep_approx),'marker','o','color','b')
legend('Jaccard distance after GA','Jaccard distance after GD')
title('Compare approximate Jaccard distance of 10 polyNbr(30-300) with 8 repeated simulation')
xlabel('Polymerase number (x3`0)')
ylabel('Approximate Jaccard distance')
ylim([0.5,1])

% ----- Jaccard distance for different amplitudes------ 
% info stored in <jaccard_dist_GD_rep> and <jaccard_dist_GA_rep_approx>
figure()
boxplot(jaccard_dist_GA_rep,'color','r')
hold on
boxplot(jaccard_dist_GD_rep,'color','b')
plot(mean(jaccard_dist_GA_rep),'marker','o','color','r')
plot(mean(jaccard_dist_GD_rep),'marker','o','color','b')
legend('Jaccard distance after GA','Jaccard distance after GD')
title('Compare Jaccard distance of 10 amplitudes(0-2.7) with 8 repeated simulation(same art_data)')
xlabel('Amplitude of noise(x0.3)')
ylabel('Jaccard distance')

figure()
boxplot(jaccard_dist_GA_rep_approx,'color','r')
hold on
boxplot(jaccard_dist_GD_rep_approx,'color','b')
plot(mean(jaccard_dist_GA_rep_approx),'marker','o','color','r')
plot(mean(jaccard_dist_GD_rep_approx),'marker','o','color','b')
legend('Jaccard distance after GA','Jaccard distance after GD')
title('Compare approximate Jaccard distance of 10 amplitudes(0-2.7) with 8 repeated simulation(same art_data)')
xlabel('Amplitude of noise(x0.3)')
ylabel('Approximate Jaccard distance')


%% ------------ Compare positions --------------
% Every row represents the positions of simulation after GD. Different rows
% have different amplitude. Blue vertical lines represent the positions of
% artificial data. If a simu position(circles) is more closed to the
% art_position, it means that our algo is more precise. Red triangles
% represent the positions where simu_position=art_postion(we found the right answer)

figure()
hold on
x = find(Pattern_poly_art_diff_amp==1);
y = 1:0.1:9;
[X,Y] = meshgrid(x,y); 
line(X,Y,'color',[0,0.7,0.9])
legend('Artificial polymerase position')
for ii = 1:length(poly_position_diff_amp_GD_rep)
    amp_i = poly_position_diff_amp_GD_rep{ii};
    p_posi = amp_i{1};
    same_position = intersect(p_posi,find(Pattern_poly_art_diff_amp==1));
    plot(p_posi,ii,'Marker','o','color',[0.01*ii,0.1*ii,0.01*ii])
    plot(same_position,ii,'Marker','^','color','red')
end
title('Polymerase position of 10 amplitudes(0-2.7) with 8 repeated simulation(same art-data)')
xlabel('Polymerase position after GD')
ylabel('Amplitude of noise(x0.3)')

%% -----------kstest--------
% ------------ KStest for different Poly Num--------------
ks_test_ga = zeros(8,10);
ks_test_gd = zeros(8,10);
for ii=1:80
        dx_a = diff(sort(art_data_diff_polyNbr(ceil(ii/8)).trans_posi_art));
        dx_s_ga = diff(sort(simu_data_diff_polyNbr_rep(ii).trans_posi));
        dx_s_gd = diff(sort(simu_data_diff_polyNbr_rep(ii).trans_posi_GDy));
        [h_ga,p_ga] = kstest2(dx_a,dx_s_ga);
        [h_ga,p_gd] = kstest2(dx_a,dx_s_gd);
        ks_test_ga(ii) = p_ga;
        ks_test_gd(ii) = p_gd;
end
figure(16)
hold on
boxplot(ks_test_ga,'color','r','Labels',30:30:300)
boxplot(ks_test_gd,'color','b','Labels',30:30:300)
plot(mean(ks_test_ga),'marker','o','color','r')
plot(mean(ks_test_gd),'marker','o','color','b')
% set(gca,'xtick',[]);
legend('KStest after GA','KStest after GD')
title('Compare KStest of 10 polyNbr(20-200) with 8 repeat)')
xlabel('polyNbr')
ylabel('p-value')


% ------------ KStest for different amplitude--------------
ks_test_ga = zeros(8,10);
ks_test_gd = zeros(8,10);
dx_a = diff(sort(art_data_diff_amp.trans_posi_art));
for ii=1:80
              
        dx_s_ga = diff(sort(simu_data_diff_amp_rep(ii).trans_posi));
        dx_s_gd = diff(sort(simu_data_diff_amp_rep(ii).trans_posi_GDy));
        [h_ga,p_ga] = kstest2(dx_a,dx_s_ga);
        [h_gd,p_gd] = kstest2(dx_a,dx_s_gd);
        ks_test_ga(ii) = p_ga;
        ks_test_gd(ii) ; p_gd;
    
end
figure()
hold on
boxplot(ks_test_ga,'color','r','Labels',[0:9]*0.3)
boxplot(ks_test_gd,'color','b','Labels',[0:9]*0.3)
plot(mean(ks_test_ga),'marker','o','color','r')
plot(mean(ks_test_gd),'marker','o','color','b')
legend('KStest after GA','KStest after GD')
title('Compare KStest of 10 amp(0-2.7) with 8 repeat)')
xlabel('amplitude of noise')
ylabel('p-value')











