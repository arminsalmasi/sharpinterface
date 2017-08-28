clc, close all, clear variables 

% % save('data.mat', 'grid1_save', ...
% %                  'grid2_save', ...
% %                  'u2_save', ...
% %                  'u1_save', ...
% %                  'reg1_save', ...
% %                  'reg2_save' , ...
% %                  'save_times',...
% %                  'save_idx',...
% %                  't' , ...
% %                  'dt', ...
% %                  'save_flag' )

d = load('data.mat');

num_plot = 5;

figure
hold on
k = 1
cnt=1
while k< d.save_idx 
          grid1_plot(:) = d.grid1_save(k,:);
          grid2_plot(:) = d.grid2_save(k,:);
          u1_plot(:,:) = d.u1_save(k,:,:);     
          u2_plot(:,:) = d.u2_save(k,:,:);     
          plot( [grid1_plot,grid2_plot], [u1_plot(1,:),u2_plot(1,:)] );
          legendCell(cnt) = cellstr(num2str(d.save_times(k)','time(s) = %0.f' ));
          k = k + floor(d.save_idx/ num_plot);
          cnt =cnt+1;
end
if k>= d.save_idx
   k = k - d.save_idx;
end    
xlim([d.reg1_save(k,1),d.reg2_save(k,2)]);
title('u-fraction of carbon - distance');
ylabel('Carbon u-fraction');
xlabel('distance meter');
legend(legendCell(:))
hold off

figure
hold on
    plot(d.save_times ,d.velocity_save);
    title('velocity of interface - time');
    xlabel('time seconds');
    xlim([0,d.t]);
    ylabel('velocity of interface m/s')   
hold off

figure
hold on
    plot(d.save_times ,d.velocity_save);
    set(gca,'Xscale', 'log')
    title('velocity of interface - time')
    xlim([0,d.t]);
    xlabel('time seconds')
    ylabel('velocity of interface m/s')   
hold off

figure
hold on
    plot(d.save_times ,d.reg1_save(:,2));
    set(gca,'Xscale', 'log')
    title('position of interface - time')
    xlim([0,d.t]);
    xlabel('time seconds')
    ylabel('position of interface meter')   
hold off




 

