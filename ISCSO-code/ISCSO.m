%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Improved Sand Cat Swarm Optimization Algorithm (ISCSO)                                      %
%                                                                                              %
%  Source codes demo version 1.0                                                               %                                                                
%                                                                                              %                                                                                                    
%  The 11th Gen Intel(R) Core(TM) i7-11700 processor with the primary frequency of 2.50GHz,    %
%  16GB memory, and the operating system of 64-bit windows 11 using matlab2023a.               %                                                       
%                                                                                              %               
%  Author and programmer: Heming Jia, Jinrui Zhang, Honghua Rao, Laith Abualigah               %                                                                      
%         E-Mail: jiaheminglucky99@126.com; Ruriuiz@gmail.com                                  %                                                   
%                                                                                              %                                                                                                  
%______________________________________________________________________________________________%
% You can simply define your cost function in a seperate file and load its handle to fobj      %
% The initial parameters that you need are:                                                    %
%__________________________________________                                                    %
% fobj = @YourCostFunction                                                                     %
% dim = number of your variables                                                               %
% T = maximum number of iterations                                                             %
% N = number of search agents                                                                  %
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n                              %
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n                              %
% If all the variables have equal lower bound you can just                                     %
% define lb and ub as two single numbers                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Best_Score,BestFit,Convergence_curve]=ISCSO(N,T,lb,ub,dim,fobj)
%% Define Parameters
BestFit=zeros(1,dim);
Best_Score=inf;
Positions=initialization(N,dim,ub,lb);
Convergence_curve=zeros(1,T);
counter=0; %Initializes the accumulator
t=0;
p=[1:360];
C = 0.01;
beta = 2.77;
while t<T
    for i=1:size(Positions,1)
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        fitness=fobj(Positions(i,:));
        if fitness<Best_Score
            Best_Score=fitness;
            BestFit=Positions(i,:);
        end
    end
    S=2; % S is maximum Sensitivity range
    rg=S-((S)*t/(T));% guides R
    for i=1:size(Positions,1)
        r=rand*rg;
        R=((2*rg)*rand)-rg; % Controls to transtion phases
        for j=1:size(Positions,2)
            %% Spiral contraction walking strategy
            teta=RouletteWheelSelection_ISCSO(p);   % Roulette options
            if((-1<=R)&&(R<=1))    % R value is between -1 and 1
                Rand_position=abs(rand*BestFit(j)-Positions(i,j));
                q = (BestFit(j)-Positions(i,j));
                e = cos(pi/2*rand)*C*(q-((q)*j/(size(Positions,2)))-1);
                Positions(i,j)=BestFit(j)+e*Rand_position*cos(teta);
                Positions(i,j)=simplebounds(Positions(i,j),lb,ub);
            else
                %% Low-frequency noise search strategy
                cp=floor(N*rand()+1);
                CandidatePosition =Positions(cp,:);
                if j>1
                    F1 = 2/(Positions(i,j)-BestFit(j-1));
                    F2 = 2/(Positions(i,j)-BestFit(j));
                    jr2 = randi(size(Positions,2));
                    F3 = 2/(Positions(i,j)-BestFit(jr2));
                    Fa = (F1+F2+F3)/3;
                    if F2>Fa
                        Positions(i,j)=F2^2*beta*r*(CandidatePosition(j)-rand*Positions(i,j));
                    else        
                         Positions(i,j)=Fa^2*beta*r*(CandidatePosition(j)-rand*Positions(i,j));
                    end
                else
                    Positions(i,j)=(0.5+rand*beta)*rg*(CandidatePosition(j)-rand*Positions(i,j));
                end
                Positions(i,j)=simplebounds(Positions(i,j),lb,ub);
            end
        end
        %% Random opposition-based learning strategy
        Xnew = lb+ub-rand*Positions(i,:);
        if fobj(Xnew)<fobj(Positions(i,:))
            Positions(i,:)=Xnew;
        end
        %% Restart strategy
        fv4=fobj(Positions(i,:));
        if fv4 < Best_Score 
            counter = 0;
        else
            counter= counter + 1;
        end
        if counter > log(t)  % Restart mechanism judgment, considering the number of iterations 1,2,5,10,15,... it/10,log(it)
            t1 = zeros(1,dim); t2 = zeros(1,dim);
            t1(1,:) = (ub-lb)*rand+lb;
            t2(1,:) = (ub+lb)*rand-Positions(i,:);
            Flag4ub=t2(1,:)>ub;
            Flag4lb=t2(1,:)<lb;
            t2(1,:)=(t2(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            if fobj(t1) < fobj(t2)
                Positions(i,:) = t1(1,:);
            else
                Positions(i,:) = t2(1,:);
            end
            counter = 0;
        end
    end
    if mod(t,50)==0  %Print the best universe details after every 50 iterations
        display(['ISCSO At iteration ', num2str(t), ' the best fitness is ', num2str(Best_Score)]);
    end
    t=t+1;
    Convergence_curve(t)=Best_Score;
end
end

%% Boundary processing
function s = simplebounds(s,lb,ub)
% Apply the lower bound
ns_tmp = s;
I = ns_tmp<lb;
ns_tmp(I) = lb;
% Apply the upper bounds
J = ns_tmp>ub;
ns_tmp(J) = ub;
% Update this new move
s = ns_tmp;
end