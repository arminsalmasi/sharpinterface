clc; clear variables, close all

format long

%% time
t = 10*3600; %simulation time
dt = 0.001;  %timestep

%inerval between saves
save_time = 100; % save each n'th second
save_flag = floor(save_time/dt); % counter on timesteps  

%% grid and domain setting
reg1 = [0, 1e-9]; % geometry of region 1 
reg2 = [reg1(2), 2e-3]; % % geometry of region 2 
np_regs = [10, 100]; % number of points in the domain [region 1, region 2]


%% composition at Boundaries
xBC = [ 7.91841E-04, 9.99208E-01  ...  
      ; 7.91841E-04, 9.99208E-01  ...
      ; 2.83247E-02, 9.71675E-01  ...
      ; 9.23190E-03, 9.90770E-01 ];

%% time loop

aclock = clock;
fprintf('%s %d-%d-%d %d:%d:%d \n', 'simulation start time:',...
    aclock(1), aclock(2), aclock(3), aclock(4), aclock(5), ...
    round(aclock(6)));
fprintf('******************* \n \n');

save_idx = 1;

%% Initiate concentration
[u1, u2, grid1, grid2, uBC_C] = do_init(reg1, reg2, np_regs, xBC);

%% pre-allocation / speed up
n_saves = floor(t / (dt * save_flag)) + 1;
velocity_save = zeros(n_saves);
grid1_save = zeros(n_saves, size(grid1,2)); 
grid2_save = zeros(n_saves, size(grid2,2)); 
u2_save = zeros(n_saves, 2, size(grid2,2));
u1_save = zeros(n_saves, 2, size(grid1,2)); 
reg1_save = zeros(n_saves, 2);
reg2_save = zeros(n_saves, 2);
save_times = zeros(n_saves);

% delete(gcp('nocreate'));
% parpool('local',2);
% parfor tstp = 1 : floor(t/dt);   

for tstp = 1 : floor(t/dt)   
    tic
    if tstp == 1
      %% save initial composition in tstp 1
      velocity=0;
      grid1_save(save_idx,:) = grid1;
      grid2_save(save_idx,:) = grid2;
      u2_save(save_idx,:,:) = u2;
      u1_save(save_idx,:,:) = u1;
      reg1_save(save_idx,:) = reg1;
      reg2_save(save_idx,:) = reg2;
      save_times(save_idx) = tstp*dt;
      save_idx = save_idx + 1;

      %% time loop info print
      fprintf('initialization loop \n');
      printloop([],tstp, dt , reg2(1), velocity , 1);
      
    else        
        %% solve PDE, calculate flux
        D2_c = 2.5e-10 * u2(1,:) + 5.9e-13; % D reg=2 el= C
        dudx2_C = diff(u2(1,:)) ./ diff(grid2); % du/dx reg=2 el= C
        flux_C_2 = -1 * D2_c(1,1:end-1) .* dudx2_C(1,1:end);
     
        %% solve flux equation, calculate interface velocity & displacement
        velocity = flux_C_2(1) / (uBC_C(3) - uBC_C(2));
        dx = velocity * dt;
      
        %% update u-fractions by flux
        h = diff(grid2);
        u2(1, 2:end-1) = u2(1, 2:end-1) + ...
                         (flux_C_2(1, 1:end-1) - flux_C_2(1, 2:end)) *...
                         dt ./ h(1,1:end-1); 
        
        %% apply boundary condition at the left side of regin 2
        u2(1, end) = u2(1, end) + flux_C_2(1, end) * dt / h(end);
      
        %% update grid
        if reg1(2)+dx>reg2(2) 
            er1 = true;
        else
            %% update regions
            er1=false;
            reg1(2) = reg1(2)+dx;
            reg2(1) = reg2(1)+dx;
            grid1(:) = linspace(reg1(1), reg1(2), np_regs(1)); 
            grid2(:) = linspace(reg2(1), reg2(2), np_regs(2)); 
        end
        
        %% save variables in keepers
        if tstp == 2
            %% loop ifo print
            printloop(toc/2, tstp, dt, reg2(1), velocity,  save_flag-tstp);
            velocity_save(1) = velocity;
        end
         if mod(tstp,save_flag) == 0 
            %% loop info print
            printloop(toc/save_flag, tstp, dt, reg2(1), velocity, ...
                     save_flag);
            %% save variables in keepers
            velocity_save(save_idx) = velocity;
            grid1_save(save_idx,:) = grid1;
            grid2_save(save_idx,:) = grid2;
            u2_save(save_idx,:,:) = u2;
            u1_save(save_idx,:,:) = u1;
            reg1_save(save_idx,:) = reg1;
            reg2_save(save_idx,:) = reg2;
            save_times(save_idx) = tstp*dt;
            save_idx = save_idx + 1;

         end  
    end
end

%% end time of simulation
bclock = clock;
fprintf('%s %d-%d-%d %d:%d:%d \n', 'simulation end time:',...
  bclock(1), bclock(2), bclock(3), bclock(4), bclock(5), round(bclock(6)));

%delete(gcp('nocreate'));

%% save results - finalize
save('data.mat', 'velocity_save', ...
                 'grid1_save', ...
                 'grid2_save', ...
                 'u2_save', ...
                 'u1_save', ...
                 'reg1_save', ...
                 'reg2_save' , ...
                 'save_times',...
                 'save_idx',...
                 't' , ...
                 'dt', ...
                 'save_flag' );

%% error handler
if er1
    fprintf('%s', ...
        "at some point Interface moved out of doamin's right baoundry \n");
end
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
function [] = printloop(Wallt_loop,timestep,dt,pos,velocity,num_to_next)
      if not(isempty(Wallt_loop))
          fprintf('Wtime/loop = %d\n', Wallt_loop )
      end
      fprintf('timestep = %d\n', timestep );
      fprintf('time = %f\n', timestep*dt );
      fprintf('interface position = %e\n', pos );
      fprintf('velocity of interface = %e\n', velocity);
      fprintf('number of timesteps to next sampling %d\n\n', num_to_next);
      fprintf('******************* \n \n');
end

function [u1, u2, grid1, grid2, uBC_C] = do_init(reg1,reg2,np_regs, xBC )

    format long
    
    %% linear grid in both regions
    grid1 = linspace(reg1(1), reg1(2), np_regs(1)); 
    grid2 = linspace(reg2(1), reg2(2), np_regs(2)); 
    
    %% concentration at BCs, u-fraction
    uBC_C =  [ xBC(1,1)/xBC(1,2), ... %right side of reg1
               xBC(2,1)/xBC(2,2), ... %left side of reg1
               xBC(3,1)/xBC(3,2), ... %right side of reg2
               xBC(4,1)/xBC(4,2)];    %left side of reg2

    %% Initial concentration profiles
    x1(1,:) = grid1 * 0 + xBC(1,1);%% x_reg1[Phase_1:alpha,element_1:C]
    x1(2,:) = grid1 * 0 + xBC(1,2);%% x_reg1[Phase_1:alpha,element_2:Fe]
    x2(1,:) = grid2 * 0 + xBC(4,1);%% x_reg2[phase_2:gamma,elemet_1:C]
    x2(2,:) = grid2 * 0 + xBC(4,2);%% x_reg2[phase_2:gamma,elemet_2:Fe]
    u1(1,:) = x1(1,:) ./ x1(2,:);  %% u_reg1[Phase_1:alpha,element_1:C]
    u1(2,:) = x1(2,:) ./ x1(2,:);  %% u_reg1[Phase_1:alpha,element_2:Fe]
    u2(1,:) = x2(1,:) ./ x2(2,:);  %% u_reg2[phase_2:gamma,elemet_1:C]
    u2(2,:) = x2(2,:) ./ x2(2,:);  %% u_reg2[Phase_2:gamma,element_2:Fe]
    
    %% fix concentrations on boundaries, u-fraction
    % region 1
    u1(1,1) = uBC_C(1);   %right side of reg1
    u1(1,end) = uBC_C(2); %left side of reg1
    % region 2
    u2(1,1) = uBC_C(3); %right side of reg2
    u2(1,end) =uBC_C(4);%left side of reg2
end