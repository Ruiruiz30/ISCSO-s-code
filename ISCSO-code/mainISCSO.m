%  Improved Sand Cat Swarm Optimization Algorithm (ISCSO)                                      %
%                                                                                              %
%  Source codes demo version 1.0                                                               %                                                                
%                                                                                              %                                                                                                    
%  The 11th Gen Intel(R) Core(TM) i7-11700 processor with the primary frequency of 2.50GHz,    %
%  16GB memory, and the operating system of 64-bit windows 11 using matlab2023a.               %                                                       
%                                                                                              %               
%  Author and programmer: Heming Jia, Jinrui Zhang, Honghua Rao, Laith Abualigah               %                                                                      
%         E-Mail: jiaheminglucky99@126.com; Ruriuiz@gmail.com                                  %     

clear all 
close all
clc
N=100;           %Number of search agents
F_name='F5';     %Name of the test function
T=200;           %Maximum number of iterations

    
[lb,ub,dim,fobj]=Get_F(F_name); %Get details of the benchmark functions
[best_fun,best_position,cuve_f]=ISCSO(N,T,lb,ub,dim,fobj); 


figure('Position',[454   445   694   297]);
subplot(1,2,1);
func_plot(F_name);     % Function plot
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([F_name,'( x_1 , x_2 )'])
subplot(1,2,2);       % Convergence plot
semilogy(cuve_f,'LineWidth',3)
xlabel('Iteration#');
ylabel('Best fitness so far');
legend('ISCSO');



display(['The best-obtained solution by ISCSO is : ', num2str(best_position)]);  
display(['The best optimal value of the objective funciton found by ISCSO is : ', num2str(best_fun)]);  