clear
clc

gen_no=500;
popusize=200;
xover_rate=1;
mutate_rate=1; % Represent each chromosome to join mutation
mutate_times=1; % Pm=1/8=0.125
stopvalue=0;
max_var=31; %L32
var_no=9; % Variables (or Factors)

a=10;
b=10; 
c=5;
lowerbound=[-b  -b  -b  0.1 0.1 0.1 -a -a -a];
            
upperbound=[b b b c c c a a a];
tic
%-------Generate and evaluate initial population-------
for k1=1:popusize
    for k2=1:var_no
        beta=rand;%random0_1;
        popu(k1,k2)=lowerbound(k2)+beta*(upperbound(k2)-lowerbound(k2));
    end %k2
    %popu(k1,1:4)=floor(popu(k1,1:4); %elements 1-4 are integers
    %popu(k1,5:8)=floor(popu(k1,5:8)*10)/10; %elements 5-8 are XX.X
    [fitness(k1)]=optimal(popu(k1,1:var_no));
end % k1

%------------Special data--------------------------------- 

% popu(1,:)=[41.286 1.2767 38.322]; %fitness=1.1928
% popu(1,:)=[50	1.252	47.754]; % 1.1790

[fitness(1)]=optimal(popu(1,1:var_no));
%---------------------------------------------------------------

%------------Calculate optimal fitness-------------------------
for n_gen=1:gen_no
    n_gen
    %------Roulette wheel method for selection operation ------
    xfitness=fitness/sum(fitness);
    cum_prob=cumsum(xfitness);
    %-----------------------------------------------------------
    xover_no=0;
    for k=1:popusize
        if rand < xover_rate
            xover_no=xover_no+1;
            tmp1=find(cum_prob > rand);
            parent1=popu(popusize-tmp1(1)+1,1:var_no);% 因望小問題,(popusize-tmp1+1)為選用finess較小者
            tmp2=find(cum_prob > rand);
            parent2=popu(popusize-tmp2(1)+1,1:var_no); % 因望小問題,(popusize-tmp2+1)為選用finess較小者
            %--------------Crossover operation-----------------------------
            xover_point=ceil(rand*var_no);
            xover_popu(xover_no*2-1,1:xover_point-1)=parent1(1:xover_point-1);
            xover_popu(xover_no*2,1:xover_point-1)=parent2(1:xover_point-1);
            
            beta=rand;%random0_1
            xover_popu(xover_no*2-1,xover_point)=(1-beta)*parent1(xover_point)+beta*parent2(xover_point);
            xover_popu(xover_no*2,xover_point)=(1-beta)*lowerbound(xover_point)+beta*upperbound(xover_point);
            
            xover_popu(xover_no*2-1,xover_point+1:var_no)=parent2(xover_point+1:var_no);
            xover_popu(xover_no*2,xover_point+1:var_no)=parent1(xover_point+1:var_no);
        end % for if
    end % for i
    %xover_popu(:,1:4)=floor(xover_popu(:,1:4)); % elements 1-4 are integers
    %xover_popu(:,5:8)=floor(xover_popu(:,5:8)*10)/10; % elements 5-8 are XX.X
    
    %------------------Taguchi Method-------------------------
    fid = fopen('L32.txt'); % 讀入欲使用之直交表
    if fid == -1
        message='Cannot open file';
        disp(message);
    else
        myData=fscanf(fid,'%d',[max_var max_var+1]); %[max_var max_var+1]是控制輸入資料陣列型態
    end
    fclose(fid);
    % L_table表使用直交表
    L_Table = myData'; %轉置成所要的直交表陣列型態
    for ta = 1 : xover_no/2
        int_rand1 = ceil(rand*2*xover_no);
        int_rand2 = ceil(rand*2*xover_no);
        taguchiTable = zeros (max_var+5,max_var+1);
        for j = 1 : max_var+1
            for k = 1 : var_no
                if L_Table(j,k) == 1
                    taguchiTable(j,k) = xover_popu(int_rand1,k);
                elseif L_Table(j,k) == 2
                    taguchiTable(j,k) = xover_popu(int_rand2,k);
                end
                popu_tmp2(k) = taguchiTable(j,k);
            end
            taguchiTable(j,max_var+1) = optimal(popu_tmp2);
            if taguchiTable(j,max_var+1) == 0.0
               taguchiTable(j,max_var+1) = 0.00001;
            end
        end
        for j = 1 : var_no
            for k = 1 : max_var+1
                if L_Table(k,j) == 1
                   taguchiTable(max_var+2,j) = taguchiTable(max_var+2,j) + 1.0 / taguchiTable(k,max_var+1);
                elseif L_Table(k,j) == 2
                   taguchiTable(max_var+3,j) = taguchiTable(max_var+3,j) + 1.0 / taguchiTable(k,max_var+1);
                end
            end
            if taguchiTable(max_var+2,j)>=taguchiTable(max_var+3,j)
               taguchiTable(max_var+4,j)=xover_popu(int_rand1,j);
            elseif taguchiTable(max_var+2,j)<taguchiTable(max_var+3,j)
                taguchiTable(max_var+4,j)=xover_popu(int_rand2,j);
            end
            popu_tmp2(j)=taguchiTable(max_var+4,j);
        end
        taguchiTable(max_var+4,max_var+1) = optimal(popu_tmp2);
        k=1;
        taguchiTable(max_var+5,max_var+1) = taguchiTable(1,max_var+1);
        for j=1:max_var+1
            if taguchiTable(j,max_var+1) < taguchiTable(max_var+5,max_var+1)
                taguchiTable(max_var+5,max_var+1) = taguchiTable(j,max_var+1);
                k=j;
            end
        end
        taguchiTable(max_var+5,1:var_no)=taguchiTable(k,1:var_no);
        if taguchiTable(max_var+4,max_var+1)>=taguchiTable(max_var+5,max_var+1)
            ta_chro(ta,1:var_no)=taguchiTable(max_var+5,1:var_no);
            fcn_vv(ta)=taguchiTable(max_var+5,max_var+1);
        elseif taguchiTable(max_var+4,max_var+1)<taguchiTable(max_var+5,max_var+1)
            ta_chro(ta,1:var_no)=taguchiTable(max_var+4,1:var_no);
            fcn_vv(ta)=taguchiTable(max_var+4,max_var+1);
        end
   end % for ta
%------------Mutation operation----------------
    for i = 1: xover_no/2
        popu_tmp1(1:var_no)=ta_chro(i,1:var_no);
        if rand < mutate_rate
            for j3 = 1:mutate_times
                int_rand1 = ceil(rand*var_no);
                beta=rand;%random0_1
                if rand > 0.5
                    popu_tmp1(int_rand1)=popu_tmp1(int_rand1)+(upperbound(int_rand1)-popu_tmp1(int_rand1))*beta;
                else
                    popu_tmp1(int_rand1)=popu_tmp1(int_rand1)-(popu_tmp1(int_rand1)-lowerbound(int_rand1))*beta;
                end
            end
        end % for if 
        % popu_tmp1(1:4)=floor(popu_tmp1(1:4)); % elements 1-4 are integers
        % popu_tmp1(5:8)=floor(popu_tmp1(5:8)*10)/10; % elements 5-8 are XX.X

        double_tmp1=optimal(popu_tmp1);
        if double_tmp1<fcn_vv(i) % if mutation has smaller fitness, replace original fitness.
            fcn_vv(i)=double_tmp1;
            ta_chro(i,1:var_no)=popu_tmp1(1:var_no);
        end % if 
    end % for i ------------------------
    %-----integrate fitness and population of popu and Ta_chro-----
    total_popu=[]; total_fitness=[];sort_total_fitness=[];sort_total_popu=[];% reset
    total_popu=[popu      % fitness
                ta_chro]; % fcn_vv
    total_fitness=[fitness fcn_vv];
    
    [sort_total_fitness,YY]=sort(total_fitness); % sorting from smaller fitness to larger fitness
    sort_total_popu=total_popu(YY,:);
    
    %-----Generate new population------
    fitness=[];popu=[];% reset
    for gg=1:popusize
       fitness(gg)=sort_total_fitness(gg);
       popu(gg,:)=sort_total_popu(gg,:);
    end % for gg-----------------------

    %----------Print fitness to display---------
    disp('fitness='); disp(fitness(1:5));
    gen_optimum(n_gen)=fitness(1); % 紀錄每一代第一個最小者)
                       
   %------------Re-initiation----------------------
   for j1=1:popusize-1
       flag=0;
       if fitness(j1)~=fitness(j1+1)
         flag=1;
       end
       if flag==1
         break;
       end
   end % for 
   if flag==0
       for k1=2:popusize % 保留第一個最小者
           for k2=1:var_no
               beta=rand;%random0_1
               popu(k1,k2)=lowerbound(k2)+beta*(upperbound(k2)-lowerbound(k2));
           end %k2
        % popu(k1,1:4)=floor(popu(k1,1:4)); % elements 1-4 are integers
        % popu(k1,5:8)=floor(popu(k1,5:8)*10)/10; % elements 5-8 are XX.X
        [fitness(k1)]=optimal(popu(k1,1:var_no));
       end % k1
    end % if flag
        % ------Stop Condition--------------------
   if gen_optimum(n_gen)<=stopvalue
      break;
   end
   
end % for n_gen
%-------------------------------------------------------
    times=popusize+(popusize*xover_rate/2)*(4+1+1)*gen_no
save('output_gen_optimum_HTGA','popusize','xover_rate','mutate_rate','n_gen','mutate_times','gen_optimum','fitness','popu','times');
t=toc/60;
fprintf('==> Computation time is (%.2f) minutes.\n\n', t);
n_gen

% popu_t=popu(1,:);
% 
% t=[25	25	25	25	30	30	35	35	35	40	40	40	40;
%   3.5	4	4	4.5	4.5	4.5	4	5	3.5	5	4	5	5;
%   0.173	0.207	0.238	0.25	0.173	0.207	0.173	0.173	0.25	0.173	0.207	0.238	0.25;
%   2.1051	1.869	1.634	0.5328	2.7427	2.3003	2.5109	2.6405	1.9633	2.1102	2.4161	1.9147	1.5193];
% 
% x=(t(1,:)-20)/(40-20);%input x  "spindle"
% y=(t(2,:)-3.5)/(5-3.5);%input y "feed"
% u=(t(3,:)-0.173)/(0.25-0.173);%input u "depth"
% z_nor=(t(4,:)-0.5328)/(2.8113-0.5328);%result
% 
% md=[x;y;u;z_nor];
% p=md';
% a11=popu_t(1);
% a12=popu_t(2);
% 
% a21=popu_t(3);
% a22=popu_t(4);
% 
% a31=popu_t(5);
% a32=popu_t(6);
% 
% b11=popu_t(7);
% b12=popu_t(8);
% 
% b21=popu_t(9);
% b22=popu_t(10);
% 
% b31=popu_t(11);
% b32=popu_t(12);
% 
% %output=wi*(px+qy+ru+s)=z
% %1**
% p111=popu_t(13);
% q111=popu_t(14);
% r111=popu_t(15);
% s111=popu_t(16);
% 
% p222=popu_t(17);
% q222=popu_t(18);
% r222=popu_t(19);
% s222=popu_t(20);
% 
% %gauss funtion
% x=p(:,1);
% m11=exp(-(x-a11).^2/b11^2/2);
% m12=exp(-(x-a12).^2/b12^2/2);
% 
% y=p(:,2);
% m21=exp(-(y-a21).^2/b21.^2/2);%gauss funtion
% m22=exp(-(y-a22).^2/b22.^2/2);
% 
% u=p(:,3);
% m31=exp(-(u-a31).^2/b31.^2/2);%gauss funtion
% m32=exp(-(u-a32).^2/b32.^2/2);
% 
% %w1**
% w111=m11.*m21.*m31;
% w222=m12.*m22.*m32;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W=w111+w222;
% 
% W111=w111./W;
% W222=w222./W;
% 
% %後鑑部參數output
% %output=wi*(px+qy+ru+s)=z
% z111=(p111*p(:,1)+q111*p(:,2)+r111*p(:,3)+s111).*W111;
% z222=(p222*p(:,1)+q222*p(:,2)+r222*p(:,3)+s222).*W222;
% 
% Z=z111+z222;
% 
% %還原預測值
% Z_exa=(2.8113-0.5328)*Z'+0.5328;%result
% 
% %誤差均方值
% z=t(4,:);
% j=sqrt(sum((z-Z_exa).^2)/(13))