% Name - Mridul Gupta
% Roll Number - 19IM30025

clear all
close all
clc

% GA Parameters
pop_size = 60;
crossover_probability = 0.8;
mutation_probability = 0.006;
number_of_generations = 60;

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

% Generating random population
n=10;
for i = 1:pop_size
    x = randi([0, 1], [1, n]);
    while sum(x.*L) > (1-K)*D
        x = randi([0, 1], [1, n]);
    end
    pop(i,:) = x;
end

% Calculating fitness for each individual of population
for i = 1:pop_size
    fitness(i) = objfunc([pop(i, :)],y,rD,D);
end

% variable for global best solution
Gbestcost = -inf;

% finding best fitness and storing in Gbestcost and best solution in best_sol 
for i = 1:pop_size
    if fitness(i) > Gbestcost
        Gbestcost = fitness(i);
        best_sol(:) = pop(i,:);
    end
end

% Main GA code
for i = 1:number_of_generations
 
    cumulative_fitness = cumsum(fitness);
    for j = 1:pop_size
        probability(j) = cumulative_fitness(j)/cumulative_fitness(pop_size);
    end
    
    % Selecting parents for crossover
    for j = 1:pop_size
        r = rand;
        for k = 1:pop_size
            if r <= probability(k)
                parent(j,:) = pop(k,:);
                break
            end
        end
    end
    
    % One point Crossover
    for j = 1:pop_size/2
        P1 = parent(j,:);
        P2 = parent(j+5,:);
        if rand <= crossover_probability
            cut_point = round(1+rand*(10-1));
            C1 = cat(2,P1(1:cut_point),P2(cut_point+1:10));
            C2 = cat(2,P2(1:cut_point),P1(cut_point+1:10));
        else
            C1 = P1;
            C2 = P2;
        end
        children(j,:) = C1;
        children(j+5,:) = C2;
    end

    % Bitwise Mutation
    for j = 1:length(children)
        if rand <= mutation_probability
            index = round(1+rand*(n-1));
            if children(j,index) == 1
                children(j,index) = 0;
            else
                children(j,index) = 1;
            end
        end
    end
    
    % Validation
    for j = 1:pop_size
        if sum(children(j).*L) <= (1-K)*D
            pop(j) = children(j);
        else
            x = randi([0, 1], [1, n]);
            while sum(x.*L) > (1-K)*D
                x = randi([0, 1], [1, n]);
            end
            pop(j,:) = x;
        end
    end
    
    % Calculate fitness of all individuals
    for j = 1:pop_size
        fitness(j) = objfunc([pop(j, :)],y,rD,D);
    end

    for j = 1:pop_size
        if fitness(j) > Gbestcost
            Gbestcost = fitness(j);
            best_sol(:) = pop(j,:);
        end
    end

    a(i) = i;
    b(i) = Gbestcost;

end

best_sol
Gbestcost
plot(a,b)
xlabel('Iteration');
ylabel('Fitness function value');