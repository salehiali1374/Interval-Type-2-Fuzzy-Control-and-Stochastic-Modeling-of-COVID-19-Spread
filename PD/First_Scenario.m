while size(dir(['*.mat']),1)<10
    clc
    clear
    close all
    warning off

    for iteration=1:10
        if exist('first_scenario_optimization_'+string(iteration)+'q5_ka54.mat')
            continue
        end
        close all

        fun = @ps;
        options.PopulationSize = 10;
        options.MaxGenerations = 500;

        [P,Fval,exitFlag,Output] = ga(fun,6,[],[],[],[],[],[],[],options);
        fprintf('The number of generations is: %d\n', Output.generations);
        fprintf('The number of function evaluations is: %d\n', Output.funccount);
        fprintf('The best function value found is: %g\n', Fval);

        %% test
        I_un_std = [];
        x_std = [];
        I_std = [];
        x_std = [];
        y_std = [];

        x = 1:1:60;
        y = 0*ones(1,60);

        % parameters
        beta = 1*ones(1,180);
        epsi = 0*ones(1,180);
        q = 0.5*ones(1,180);
        delta = 1*ones(1,180);
        kappa = 0.54*ones(1,180);
        rho = 0.1*ones(1,180);
        z = 0.02*ones(1,180);
        eta = 0.3*ones(1,180);
        alpha = 0.3*ones(1,180);
        feta = 0.965*ones(1,180);

        % control signals

        vacci_rate =0;
        isola_rate =0;

        % initialization
        S(1) = 0.9;
        E(1) = 0;
        I(1) = 0.1;
        A(1) = 0;
        R(1) = 0;

        N(1) = S(1) + E(1) + I(1) + A(1) + R(1);


        %% uncontrolled

        % parameters
        beta_un = 1*ones(1,180);
        eps_un = 0*ones(1,180);
        q_un = 0.5*ones(1,180);
        delta_un = 1*ones(1,180);
        kappa_un = 0.54*ones(1,180);
        rho_un = 0.1*ones(1,180);
        z_un = 0.02*ones(1,180);
        eta_un = 0.3*ones(1,180);
        alpha_un = 0.3*ones(1,180);
        feta_un = 0.965*ones(1,180);


        % control signals
        vacci_rate_un = 0;
        isola_rate_un = 0;

        % initialization

        S_un(1) = 0.9;
        E_un(1) = 0;
        I_un(1) = 0.1;
        A_un(1) = 0;
        R_un(1) = 0;

        N_un(1) = S_un(1) + E_un(1) + I_un(1) + A_un(1) + R_un(1);



        for i=2:60
            epsS = 0.001*S(i-1)*randn;
            epsE = 0.001*E(i-1)*randn;
            epsI = 0.001*I(i-1)*randn;
            epsA = 0.001*A(i-1)*randn;
            epsR = 0.001*R(i-1)*randn;


            epsbeta = max(0,0.01*beta(i-1)*randn);
            epseps= max(0,0.01*eps(i-1)*randn);
            epsq= max(0,0.01*q(i-1)*randn);
            epsdelta = max(0,0.01*delta(i-1)*randn);
            epskappa = max(0,0.01*kappa(i-1)*randn);
            epsrho = max(0,0.01*rho(i-1)*randn);
            epsz = max(0,0.01*z(i-1)*randn);
            epseta = max(0,0.01*eta(i-1)*randn);
            epsalpha = max(0,0.01*alpha(i-1)*randn);
            epsfeta = max(0,0.01*feta(i-1)*randn);
            epsisola_rate = max(0,0.01*isola_rate(i-1)*randn);
            epsvacci_rate = max(0,0.01*vacci_rate(i-1)*randn);

            S(i) = (S(i-1)+epsS) - (beta(i-1)+epsbeta)*((eps(i-1)+epseps)*(E(i-1)+epsE)+(1-(q(i-1)+epsq))*(I(i-1)+epsI)+(delta(i-1)+epsdelta)*(A(i-1)+epsA))...
                *(S(i-1)+epsS) - (S(i-1)+epsS)*(vacci_rate(i-1)+epsvacci_rate);
            E(i) = (E(i-1)+epsE) + (beta(i-1)+epsbeta)*((eps(i-1)+epseps)*(E(i-1)+epsE)+(1-(q(i-1)+epsq))*(I(i-1)+epsI)+(delta(i-1)+epsdelta)*(A(i-1)+epsA))...
                *(S(i-1)+epsS)-(kappa(i-1)+epskappa)*(E(i-1)+epsE);
            I(i) = (I(i-1)+epsI) + (rho(i-1)+epsrho)*(kappa(i-1)+epskappa)*(E(i-1)+epsE)+(1-(z(i-1)-epsz))*(eta(i-1)+epseta)*(A(i-1)+epsA)...
                -(alpha(i-1)+epsalpha)*(I(i-1)+epsI) - (isola_rate(i-1)+epsisola_rate)*(I(i-1)+epsI);
            A(i) = (A(i-1)+epsA) + (1-(rho(i-1)+epsrho))*(kappa(i-1)+epskappa)*(E(i-1)+epsE)-(eta(i-1)+epseta)*(A(i-1)+epsA);
            R(i) = (R(i-1)+epsR) + (z(i-1)+epsz)*(eta(i-1)+epseta)*(A(i-1)+epsA)+(feta(i-1)+epsfeta)*(alpha(i-1)+epsalpha)*(I(i-1)+epsI);

            n(i) = epsS+epsE+epsI+epsA+epsR;
            N(i) = S(i)+E(i)+I(i)+A(i)+R(i)+n(i);

            input_vector = [I(i) - y(i) ((I(i) - y(i))-(I(i-1) - y(i-1))) ((I(i) - y(i))+(I(i-1) - y(i-1)))/2]/N(i);
            input_vector(abs(input_vector)<0.01) = 0;
            evalpd_result = evalpd(P,input_vector);
            vacci_rate(i) = max(0,evalpd_result(1));
            isola_rate(i) = max(0,evalpd_result(2));

            vacci_rate(i) = min(1,vacci_rate(i));
            isola_rate(i) = min(1,isola_rate(i));

            %% uncontrolled
            epsS_un = 0.001*S_un(i-1)*randn;
            epsE_un = 0.001*E_un(i-1)*randn;
            epsI_un = 0.001*I_un(i-1)*randn;
            epsA_un = 0.001*A_un(i-1)*randn;
            epsR_un = 0.001*R_un(i-1)*randn;

            epsbeta_un = max(0,0.01*beta_un(i-1)*randn);
            epseps_un= max(0,0.01*eps_un(i-1)*randn);
            epsq_un= max(0,0.01*q_un(i-1)*randn);
            epsdelta_un = max(0,0.01*delta_un(i-1)*randn);
            epskappa_un = max(0,0.01*kappa_un(i-1)*randn);
            epsrho_un = max(0,0.01*rho_un(i-1)*randn);
            epsz_un = max(0,0.01*z_un(i-1)*randn);
            epseta_un = max(0,0.01*eta_un(i-1)*randn);
            epsalpha_un = max(0,0.01*alpha_un(i-1)*randn);
            epsfeta_un = max(0,0.01*feta_un(i-1)*randn);
            epsisola_rate_un = max(0,0.01*isola_rate_un(i-1)*randn);
            epsvacci_rate_un = max(0,0.01*vacci_rate_un(i-1)*randn);

            S_un(i) = (S_un(i-1)+epsS_un) - (beta_un(i-1)+epsbeta_un)*((eps_un(i-1)+epseps_un)*(E_un(i-1)+epsE_un)+(1-(q_un(i-1)+epsq_un))*(I_un(i-1)+epsI_un)+(delta_un(i-1)+epsdelta_un)*(A_un(i-1)+epsA_un))...
                *(S_un(i-1)+epsS_un) - (S_un(i-1)+epsS_un)*(vacci_rate_un(i-1)+epsvacci_rate_un);
            E_un(i) = (E_un(i-1)+epsE_un) + (beta_un(i-1)+epsbeta_un)*((eps_un(i-1)+epseps_un)*(E_un(i-1)+epsE_un)+(1-(q_un(i-1)+epsq_un))*(I_un(i-1)+epsI_un)+(delta_un(i-1)+epsdelta_un)*(A_un(i-1)+epsA_un))...
                *(S_un(i-1)+epsS_un)-(kappa_un(i-1)+epskappa_un)*(E_un(i-1)+epsE_un);
            I_un(i) = (I_un(i-1)+epsI_un) + (rho_un(i-1)+epsrho_un)*(kappa_un(i-1)+epskappa_un)*(E_un(i-1)+epsE_un)+(1-(z_un(i-1)-epsz_un))*(eta_un(i-1)+epseta_un)*(A_un(i-1)+epsA_un)...
                -(alpha_un(i-1)+epsalpha_un)*(I_un(i-1)+epsI_un) - (isola_rate_un(i-1)+epsisola_rate_un)*(I_un(i-1)+epsI_un);
            A_un(i) = (A_un(i-1)+epsA_un) + (1-(rho_un(i-1)+epsrho_un))*(kappa_un(i-1)+epskappa_un)*(E_un(i-1)+epsE_un)-(eta_un(i-1)+epseta_un)*(A_un(i-1)+epsA_un);
            R_un(i) = (R_un(i-1)+epsR_un) + (z_un(i-1)+epsz_un)*(eta_un(i-1)+epseta_un)*(A_un(i-1)+epsA_un)+(feta_un(i-1)+epsfeta_un)*(alpha_un(i-1)+epsalpha_un)*(I_un(i-1)+epsI_un);
            R_un(i)
            n_un(i) = epsS_un+epsE_un+epsI_un+epsA_un+epsR_un;
            N_un(i) = S_un(i)+E_un(i)+I_un(i)+A_un(i)+R_un(i)+n_un(i);

            vacci_rate_un(i) = 0;
            isola_rate_un(i) = 0;
        end

        I_un_std = [I_un_std ; I_un];
        x_std = [x_std; x];
        I_std = [I_std; I];
        x_std = [x_std; x];
        y_std = [y_std; y];

        plot(x,I_un,x,I,x,y)

        hold on
        plot(vacci_rate)
        plot(isola_rate)
        legend('Uncontrolled I','Controlled I','Desired', 'Vaccination', 'Isolation')

        mse_storage(iteration) = mse(y-I)
        if mse_storage(iteration)~=Inf & mse_storage(iteration)<0.5
            save('first_scenario_optimization_'+string(iteration)+'q5_ka54.mat')
        end
    end

    %%
    % train
end
function cost = ps(P)
x = 1:1:60;
y = 0*ones(1,60);

% parameters
beta = 1*ones(1,180);
epsi = 0*ones(1,180);
q = 0.5*ones(1,180);
delta = 1*ones(1,180);
kappa = 0.54*ones(1,180);
rho = 0.1*ones(1,180);
z = 0.02*ones(1,180);
eta = 0.3*ones(1,180);
alpha = 0.3*ones(1,180);
feta = 0.965*ones(1,180);

% control signals

vacci_rate = 0;
isola_rate = 0;

% initialization
S(1) = 0.9;
E(1) = 0;
I(1) = 0.1;
A(1) = 0;
R(1) = 0;

N(1) = S(1) + E(1) + I(1) + A(1) + R(1);


for i=2:60
    epsS = 0.001*S(i-1)*randn;
    epsE = 0.001*E(i-1)*randn;
    epsI = 0.001*I(i-1)*randn;
    epsA = 0.001*A(i-1)*randn;
    epsR = 0.001*R(i-1)*randn;


    epsbeta = max(0,0.01*beta(i-1)*randn);
    epseps= max(0,0.01*eps(i-1)*randn);
    epsq= max(0,0.01*q(i-1)*randn);
    epsdelta = max(0,0.01*delta(i-1)*randn);
    epskappa = max(0,0.01*kappa(i-1)*randn);
    epsrho = max(0,0.01*rho(i-1)*randn);
    epsz = max(0,0.01*z(i-1)*randn);
    epseta = max(0,0.01*eta(i-1)*randn);
    epsalpha = max(0,0.01*alpha(i-1)*randn);
    epsfeta = max(0,0.01*feta(i-1)*randn);
    epsisola_rate = max(0,0.01*isola_rate(i-1)*randn);
    epsvacci_rate = max(0,0.01*vacci_rate(i-1)*randn);


    S(i) = (S(i-1)+epsS) - (beta(i-1)+epsbeta)*((eps(i-1)+epseps)*(E(i-1)+epsE)+(1-(q(i-1)+epsq))*(I(i-1)+epsI)+(delta(i-1)+epsdelta)*(A(i-1)+epsA))...
        *(S(i-1)+epsS) - (S(i-1)+epsS)*(vacci_rate(i-1)+epsvacci_rate);
    E(i) = (E(i-1)+epsE) + (beta(i-1)+epsbeta)*((eps(i-1)+epseps)*(E(i-1)+epsE)+(1-(q(i-1)+epsq))*(I(i-1)+epsI)+(delta(i-1)+epsdelta)*(A(i-1)+epsA))...
        *(S(i-1)+epsS)-(kappa(i-1)+epskappa)*(E(i-1)+epsE);
    I(i) = (I(i-1)+epsI) + (rho(i-1)+epsrho)*(kappa(i-1)+epskappa)*(E(i-1)+epsE)+(1-(z(i-1)-epsz))*(eta(i-1)+epseta)*(A(i-1)+epsA)...
        -(alpha(i-1)+epsalpha)*(I(i-1)+epsI) - (isola_rate(i-1)+epsisola_rate)*(I(i-1)+epsI);
    A(i) = (A(i-1)+epsA) + (1-(rho(i-1)+epsrho))*(kappa(i-1)+epskappa)*(E(i-1)+epsE)-(eta(i-1)+epseta)*(A(i-1)+epsA);
    R(i) = (R(i-1)+epsR) + (z(i-1)+epsz)*(eta(i-1)+epseta)*(A(i-1)+epsA)+(feta(i-1)+epsfeta)*(alpha(i-1)+epsalpha)*(I(i-1)+epsI);

    n(i) = epsS+epsE+epsI+epsA+epsR;
    N(i) = S(i)+E(i)+I(i)+A(i)+R(i)+n(i);

    input_vector = [I(i) - y(i) ((I(i) - y(i))-(I(i-1) - y(i-1))) ((I(i) - y(i))+(I(i-1) - y(i-1)))/2]/N(i);
    input_vector(abs(input_vector)<0.01) = 0;
    evalpd_result = evalpd(P,input_vector);
    vacci_rate(i) = max(0,evalpd_result(1));
    isola_rate(i) = max(0,evalpd_result(2));

    vacci_rate(i) = min(1,vacci_rate(i));
    isola_rate(i) = min(1,isola_rate(i));

end

cost = 10*sqrt(sum((I-y).^2)) + sum(isola_rate((I-y)>0.01)) + sum(vacci_rate((I-y)>0.0001)) + sum(I(y==0))+sum((I-0.1)>0)+sum(vacci_rate)+sum(isola_rate);% + sum(abs(y-I)<0.0001); % cost function

end
function evalpd_result = evalpd(P,input_vector)
evalpd_result(1) = sum(P(1:3).*input_vector);
evalpd_result(2) = sum(P(4:6).*input_vector);
end