% Script to stochastic SIR model with monte carlo method
clear
close all
randn('state',3);
rand('state',3);

t0 = 0;
dt = 0.01;
N = 500;

% Number of Monte Carlo samples.
nsamples = 2000;

%Generate the samples of the random coefficient.
y_mc(1,:) = 0.7+0.1*rand(1,nsamples);
y_mc(2,:) = 0.1+0.1*rand(1,nsamples);
y_mc(3,:) = 1-(y_mc(1,:)+y_mc(2,:));

alphs_mc = 2 + 4*rand(1,nsamples);
Y_mc = zeros(3,N+1,nsamples);

% For each sample, run a deterministic simulation up to time T.
for n=1:nsamples
    % Run the simulation up to time T.
    yprime_mc = @(tt, yy) sir_rhs(yy, alphs_mc(n));
    % Time integration and output
    Y1_mc = rk4(y_mc(:,n)', yprime_mc, t0, dt, N);

    % Store the results of the simulation
    Y_mc(:,:,n) = Y1_mc;
end

%% MCS-const statistics  
nbin = 100; % 500-0.5
figure; set(0,'defaultaxesfontsize',10);  
for j=1:N
    [nIr_co1,bin_co1]=hist(Y_mc(1,j,:),nbin);  v_co1=nIr_co1/(bin_co1(2)-bin_co1(1))/sum(nIr_co1);% MCS pdf
    [nIr_co2,bin_co2]=hist(Y_mc(2,j,:),nbin);  v_co2=nIr_co2/(bin_co2(2)-bin_co2(1))/sum(nIr_co2);% MCS pdf
    [nIr_co3,bin_co3]=hist(Y_mc(3,j,:),nbin);  v_co3=nIr_co3/(bin_co3(2)-bin_co3(1))/sum(nIr_co3);% MCS pdf
    plot(bin_co1, v_co1, 'r-', bin_co2, v_co2, 'g-',bin_co3, v_co3, 'b-');
    
    set(gca,'FontSize',20);
    set(xlabel('Ratio', 'Fontsize', 25), 'interpreter', 'latex');
    set(ylabel('PDF', 'Fontsize', 25), 'interpreter', 'latex','Rotation', 90);
    h=legend('$S$', '$I$', '$R$');legend('boxoff');set(h, 'interpreter', 'latex','FontSize',20);
    pause(0.1);
end

