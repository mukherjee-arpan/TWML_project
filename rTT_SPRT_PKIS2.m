clc
clear

%% read data: uncomment to extract 4 kinase inhibitors from the dataset. 
%  The dataset has been named "drug_data.csv", and can be downloaded at:
%  https://pubmed.ncbi.nlm.nih.gov/28767711/


% data = csvread('drug_data.csv',1,0);
% data = data(data~=0);
% data = data/100;
% data = log(1 - data);
% K = 4;
% rem = data;
% mu = zeros(1,K);
% s = 1;
% while s<=K
%     samp = randsample(rem,1);
%     if min(abs(max(mu)-samp))<0.12
%         continue
%     else
%         rem = setdiff(rem,samp);
%         mu(s) = samp;
%         s = s+1;
%     end
% end
% mu = sort(mu,'descend');

%% The sampled data used for current simulations
mu = [-0.1278   -0.2877   -0.9676   -1.4697];
delta = [0.05, 0.2];
rSPRT = zeros(1,length(delta));
sigma = 1;
%delta = 0.01;
alpha_n = 1;
beta_n = 1;
eps = 0.1;
K = length(mu);
%M = 1000;
avg = 200;
for delta_iter = 1:length(delta)
%% rTT-SPRT
tau_SPRT = zeros(1,avg);
N_SPRT = zeros(1,K);
beta_SPRT = 0.5;
count = 0;
list_count = {};
time = zeros(1,avg);
correct = zeros(1,avg);
reward_list = cell(K,1);
for arm = 1:K
    reward_list{arm} = [];
end
for avg_iter = 1:avg
    tic
    % Counter of arm-pulls
    N = zeros(1,K);
    % Initialization: 
    mu_n = 0.01*ones(1,K);
    t = 1;
    LLR = zeros(1,K);
    J_val = 0;
 
    while (J_val<=log(2*t*(K-1)/delta(delta_iter)))
        if t <=K
            idx_arm = t;
            t = t+1;
            N(idx_arm) = N(idx_arm) + 1;
        elseif min(N) <= sqrt(t/K)
            [~,idx_arm] = min(N);
        else
            % Find current best
            [~,I] = max(mu_n);
            % compute GLRTs with respect to current best
            LLR(I) = inf;
            for arm = 1:K
                if arm~=I
                    mu_med = (N(arm)*mu_n(arm) + N(I)*mu_n(I))/(N(arm)+N(I));
                    %KL_I = mu_n(I)*log2(mu_n(I)/mu_med) + (1-mu_n(I))*log2((1-mu_n(I))/(1-mu_med));
                    
                    KL_I = (mu_n(I)-mu_med)^2/(2*sigma^2);
                    KL_arm = (mu_n(arm)-mu_med)^2/(2*sigma^2);

                    %KL_arm = mu_n(arm)*log2(mu_n(arm)/mu_med) + (1-mu_n(arm))*log2((1-mu_n(arm))/(1-mu_med));
                    LLR(arm) = (N(I)*KL_I + N(arm)*KL_arm);
                end
            end
            [J_val,J] = min(LLR);
 
 
 
            B = binornd(1,beta_SPRT);
            if B == 0
                %[~,idx_arm] = min(N);
                idx_arm = J;
            else
                idx_arm = I; 
            end
            t = t+1;
        end

        % Adversary flipping coin and contaminating reward
        D = binornd(1,eps);
        if D == 1
            x = 0.3+rand();
        else
            x = normrnd(mu(idx_arm),sigma);
        end
        reward_list{idx_arm} = [reward_list{idx_arm}, x];
        % Update the counter for the arm
        N(idx_arm) = N(idx_arm) + 1;
        % assign scores to each sample

        score = zeros(length(reward_list{idx_arm}),1);
        med = median(reward_list{idx_arm});
        for score_iter = 1:length(score)
            score(score_iter) = abs(reward_list{idx_arm}(score_iter)-med)/sigma;
        end

        % find probable outliers
        [~,idx_score] = sort(score,'ascend');
        good_samples = reward_list{idx_arm}(idx_score(1:ceil(N(idx_arm)*eps)));
        
        % Update posterior mean of the selected arm
        mu_n(idx_arm) = mean(good_samples);
    end
    [~,Istar] = max(mu_n);
    correct(avg_iter) = 1 == Istar;
    % Sample complexity
    tau_SPRT(avg_iter) = tau_SPRT(avg_iter) + sum(N);
    
    time(avg_iter) = toc;
    avg_iter
end
% Average Sample Complexity
rSPRT(delta_iter) = mean(tau_SPRT);
end
