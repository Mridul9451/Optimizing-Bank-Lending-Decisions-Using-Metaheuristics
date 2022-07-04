%Name - Mridul Gupta
%Roll Number - 19IM30025

clear all
close all
clc

%GA Parameters
pop_size = 60;
crossover_probability = 0.8;
mutation_probability = 0.006;

%PSO Parameters
n=10;               % no of customers
w=0.3;              % inertia weight
wdamp = 0.99;       % inertia deamping
c1=0.6;             % acceleration factor P_best
c2=1;               % acceleration factor G_best
max_iteration=60;   % maximum number of iteration in each run
maxrun=60;          % maximum number of runs

D = 60;
K = 0.15;
L = [10, 25, 4, 11, 18, 3, 17, 15, 9, 10];
rL = [0.021, 0.022, 0.021, 0.027, 0.025, 0.026, 0.023, 0.021, 0.028, 0.022];
lambda = [0.0002, 0.0058, 0.0001, 0.0003, 0.0024, 0.0002, 0.0058, 0.0002, 0.001, 0.001];
T = (1-K)*D-L;
rD = 0.009;
rT = 0.01;
for i = 1:length(L)
    y(i) = rL(i)*L(i) - lambda(i) + rT*T(i) - lambda(i);
end

%main loop
for run=1:maxrun

    for i=1:pop_size
        for j=1:n
            x0(i,j)=round(rand());
        end
    end

    %randomly generated population and check for validation
    for i = 1:pop_size     
        if sum(x0(i,:).*L) > (1-K)*D
            q = randi([0, 1], [1, n]);
            while sum(q.*L) > (1-K)*D
                q = randi([0, 1], [1, n]);
            end
            x0(i,:) = q;
        end
    end

    x = x0;                       % initial population
    v = 0.1*x0;                   % initial velocity
    
    %Calculating fitness
    for i=1:pop_size
        f0(i,1) = objfunc(x0(i,:),y,rD,D);
    end
    
    [fmax0,index0] = max(f0);     % finding highest fitness value with index
    lbest = x0;                   % initial lbest
    gbest = x0(index0,:);         % initial gbest
    no_of_iteration = 1;

    while no_of_iteration <= max_iteration

        w=w*wdamp;              % update inertia weight
        
        % calculating velocity
        for i=1:pop_size
            for j=1:n
                v(i,j)=w*v(i,j)+c1*rand()*(lbest(i,j)-x(i,j))+c2*rand()*(gbest(1,j)-x(i,j));
            end
        end
        
        % updating position
        for i=1:pop_size
            for j=1:n
                x(i,j)=round(x(i,j)+v(i,j));
            end
        end

        % Check for boundary violations
        for i=1:pop_size
            for j=1:n
                if x(i,j)<0
                    x(i,j)=0;
                elseif x(i,j)>1
                    x(i,j)=1;
                end
            end
        end

        % One point Crossover for random 10 pair of individuals
        for i = 1:n
            r1 = round(1+rand()*(pop_size-1));
            r2 = round(1+rand()*(pop_size-1));
            P1 = x(r1,:);
            P2 = x(r2,:);
            if rand <= crossover_probability
                cut_point = round(1+rand*(n-1));
                C1 = cat(2,P1(1:cut_point),P2(cut_point+1:10));
                C2 = cat(2,P2(1:cut_point),P1(cut_point+1:10));
            else
                C1 = P1;
                C2 = P2;
            end
            x(r1,:) = C1;
            x(r2,:) = C2;
        end

        % Bitwise Mutation for random 10 individuals
        for i = 1:n
            r = round(1+rand()*(pop_size-1));
            if rand <= mutation_probability
                index = round(1+rand*(n-1));
                if x(r,index) == 1
                    x(r,index) = 0;
                else
                    x(r,index) = 1;
                end
            end
        end
        
        %validation
        for i = 1:pop_size     
            if sum(x(i,:).*L) > (1-K)*D
                q = randi([0, 1], [1, n]);
                while sum(q.*L) > (1-K)*D
                    q = randi([0, 1], [1, n]);
                end
                x(i,:) = q;
            end
        end

        % evaluating fitness
        for i=1:pop_size
            f(i,1) = objfunc(x(i,:),y,rD,D);
        end

        % updating lbest and fitness
        for i=1:pop_size
            if f(i,1)>f0(i,1)
                lbest(i,:) = x(i,:);
                f0(i,1) = f(i,1);
            end
        end

        [fmax,index] = max(f0);             % finding out the best particle
       
        % updating gbest and best fitness
        if fmax>fmax0
            gbest=lbest(index,:);
            fmax0=fmax;
        end
        
        a(no_of_iteration) = no_of_iteration;
        b(no_of_iteration) = objfunc(gbest,y,rD,D);
        no_of_iteration=no_of_iteration+1;
    end
end

gbest
best_fitness_value = objfunc(gbest,y,rD,D)
plot(a,b)
xlabel('Iteration');
ylabel('Fitness function value');