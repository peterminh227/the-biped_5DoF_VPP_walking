function [] = make_plots(t,x,u,y,dy,f,step_t,doplots,doprints,print_type)


if nargin < 10
  print_type = 'ps';
end

if nargin < 9
  doprints = 0;
end

if nargin < 8
  doplots = ones(6,1);
end

if length(doplots) < 6
  error('Need to specify whether to plot EACH plot.')
end

figno = 2;
plot_num = 0;

offset = [50 50 0 0];
screen_position = [30 100 300 400];
paper_position = [.25 2.5 8 6];

plot_num = plot_num + 1;
if doplots(plot_num)
  fig_hl=figure(figno); figno = figno + 1;
  set(fig_hl,'Position', screen_position);
  screen_position = screen_position + offset;
  set(fig_hl,'PaperPosition',paper_position)

  subplot(3,1,1)
  plot(t,x(:,1),'-',t,x(:,2),'--')
  legend('q_1 (rad)','q_2 (rad)',-1)
  grid
  subplot(3,1,2)
  plot(t,x(:,3),'-',t,x(:,4),'--')
  legend('q_3 (rad)','q_4 (rad)',-1)
  grid
  subplot(3,1,3)
  plot(t,x(:,5))
  legend('q_5 (rad)',-1)
  xlabel('time (sec)')
  grid
  if doprints
    name = 'positions';
    switch print_type
     case 'ps'
      print('-dps',['plots/',name,'.ps']);
     case 'png'
      print('-dpng',['plots/',name,'.png']);
     case 'fig'
      saveas(fig_hl,['plots/',name,'.fig'],'fig');
    end
  end
end

plot_num = plot_num + 1;
if doplots(plot_num)
  fig_hl=figure(figno); figno = figno + 1;
  set(fig_hl,'Position', screen_position);
  screen_position = screen_position + offset;
  set(fig_hl,'PaperPosition',paper_position)

  subplot(3,1,1)
  plot(t,x(:,6),'-',t,x(:,7),'--')
  legend('q_{1,dot} (rad/sec)','q_{2,dot} (rad/sec)')
  grid
  axis([min(t) max(t) min(min(x(:,6:7)))*1.1 max(max(x(:,6:7)))*1.1])
  subplot(3,1,2)
  plot(t,x(:,8),'-',t,x(:,9),'--')
  legend('q_{3,dot} (rad/sec)','q_{4,dot} (rad/sec)')
  grid
  axis([min(t) max(t) min(min(x(:,8:9)))*1.1 max(max(x(:,8:9)))*1.1])
  subplot(3,1,3)
  plot(t,x(:,10))
  legend('q_{5,dot} (rad/sec)')
  grid
  xlabel('time (sec)')
  if doprints
    name = 'velocities';
    switch print_type
     case 'ps'
      print('-dps',['plots/',name,'.ps'])
     case 'png'
      print('-dpng',['plots/',name,'.png']);
     case 'fig'
      saveas(fig_hl,['plots/',name,'.fig'],'fig');
    end
  end
end

plot_num = plot_num + 1;
if doplots(plot_num)
  fig_hl=figure(figno); figno = figno + 1;
  set(fig_hl,'Position', screen_position);
  screen_position = screen_position + offset;
  set(fig_hl,'PaperPosition',paper_position)

  subplot(2,1,1)
  plot(t,u(:,1),'-',t,u(:,2),'--')
  legend('u_1 (Nm)','u_2 (Nm)')
  title('Hip controls')
  grid
  axis([min(t) max(t) min(min(u(:,1:2)))*1.1 max(max(u(:,1:2)))*1.1])
  subplot(2,1,2)
  plot(t,u(:,3),'-',t,u(:,4),'--')
  legend('u_3 (Nm)','u_4 (Nm)')
  title('Knee controls')
  xlabel('time (sec)')
  grid
  axis([min(t) max(t) min(min(u(:,3:4)))*1.1 max(max(u(:,3:4)))*1.1])
  if doprints
    name = 'control';
    switch print_type
     case 'ps'
      print('-dps',['plots/',name,'.ps'])
     case 'png'
      print('-dpng',['plots/',name,'.png']);
     case 'fig'
      saveas(fig_hl,['plots/',name,'.fig'],'fig');
    end
  end
end

plot_num = plot_num + 1;
if doplots(plot_num)
  fig_hl=figure(figno); figno = figno + 1;
  set(fig_hl,'Position', screen_position);
  screen_position = screen_position + offset;
  set(fig_hl,'PaperPosition',paper_position)

  nn = 4; mm = 1;
  
  % positions
  subplot(nn,mm,1)
  plot(t,y(:,1))
  title('Outputs (note: integrator tolerance set to 10^{-5})')
  ylabel('y_1 (rad)')
  grid
  if 0
    axis([min(t) max(t) -1e-4 1e-4]);
  else
    axis([min(t) max(t) min(min(y(:,1)))*1.1 max(max(y(:,1)))*1.1]);
  end
  subplot(nn,mm,mm+1)
  plot(t,y(:,2))
  ylabel('y_2 (rad)')
  grid
  if 0
    axis([min(t) max(t) -1e-4 1e-4]);
  else
    axis([min(t) max(t) min(min(y(:,2)))*1.1 max(max(y(:,2)))*1.1]);
  end
  subplot(nn,mm,2*mm+1)
  plot(t,y(:,3))
  ylabel('y_3 (rad)')
  grid
  if 0
    axis([min(t) max(t) -1e-4 1e-4]);
  else
    axis([min(t) max(t) min(min(y(:,3))) max(max(y(:,3)))*1.1]);
  end
  subplot(nn,mm,3*mm+1)
  plot(t,y(:,4))
  ylabel('y_4 (rad)')
  xlabel('time (sec)')
  grid
  if 0
    axis([min(t) max(t) -1e-4 1e-4]);
  else 
    axis([min(t) max(t) min(min(y(:,4))) max(max(y(:,4)))])
  end
 
  if mm == 2
    % velocities
    subplot(nn,mm,2)
    plot(t,dy(:,1))
    title('Output time derivatives')
    ylabel('dy_1 (q1-q1d)')
    xlabel('time (sec)')
    grid
    if 0
      axis([min(t) max(t) -1e-4 1e-4]);
    else
      axis([min(t) max(t) min(min(dy(:,1)))*1.1 max(max(dy(:,1)))*1.1]);
    end
    subplot(nn,mm,4)
    plot(t,dy(:,2))
    ylabel('dy_2 (d1+d2)')
    xlabel('time (sec)')
    grid
    if 0
      axis([min(t) max(t) -1e-4 1e-4]);
    else
      axis([min(t) max(t) min(min(dy(:,2)))*1.1 max(max(dy(:,2)))*1.1]);
    end
    subplot(nn,mm,6)
    plot(t,dy(:,3))
    ylabel('stance leg knee')
    xlabel('time (sec)')
    grid
    if 0
      axis([min(t) max(t) -1e-4 1e-4]);
    else
      axis([min(t) max(t) min(min(dy(:,3))) max(max(dy(:,3)))*1.1]);
    end
    subplot(nn,mm,8)
    plot(t,dy(:,4))
    ylabel('swing leg knee')
    xlabel('time (sec)')
    grid
    if 0
      axis([min(t) max(t) -1e-4 1e-4]);
    else 
      axis([min(t) max(t) min(min(dy(:,4))) max(max(dy(:,4)))])
    end
  end
  if doprints
    name = 'outputs';
    switch print_type
     case 'ps'
      print('-dps',['plots/',name,'.ps'])
     case 'png'
      print('-dpng',['plots/',name,'.png']);
     case 'fig'
      saveas(fig_hl,['plots/',name,'.fig'],'fig');
    end
  end
end

plot_num = plot_num + 1;
if doplots(plot_num)
  modelP; % get model parameters

  fig_hl=figure(figno); figno = figno + 1;
  set(fig_hl,'Position', screen_position);
  screen_position = screen_position + offset;
  set(fig_hl,'PaperPosition',paper_position)
  
  subplot(2,1,1)
  plot(t,hip_pos(x))
  ylabel('Hip horizontal position')
  xlabel('time (sec)')
  subplot(2,1,2)
  plot(t,swing_foot_height(x))
  grid
  xlabel('time (sec)')
  ylabel('Foot height (m)')
  axis([min(t) max(t) 0 max(swing_foot_height(x))*1.1])
  if doprints
    name = 'foot_hip_traj';
    switch print_type
     case 'ps'
      print('-dps',['plots/',name,'.ps'])
     case 'png'
      print('-dpng',['plots/',name,'.png']);
     case 'fig'
      saveas(fig_hl,['plots/',name,'.fig'],'fig');
    end
  end
end


plot_num = plot_num + 1;
if doplots(plot_num)
  modelP; % get model parameters

  fig_hl=figure(figno); figno = figno + 1;
  set(fig_hl,'Position', screen_position);
  screen_position = screen_position + offset;
  set(fig_hl,'PaperPosition',paper_position)
  
  [hHh,hHv]=hip_pos(x);
  [hCh,hCv] = CoM_pos(x);
  subplot(2,1,1)
  plot(t,hHv,t,hCv,'-')
  grid
  xlabel('time (sec)')
  title('Hip/CoM height (m)')
  axis([min(t),	max(t), min(hHv), max(hCv)]);
  %
  %
  subplot(2,1,2)
  plot(t,x(:,3),'-',t,x(:,4),'--')
  legend('stance leg','swing leg')
  grid
  xlabel('time (sec)')
  title('Relative knee angles (rad)')
  axis([min(t),	max(t),	min(min(x(:,3:4))), max(max(x(:,3:4)))]);
  
  if doprints
    name = 'hip_height';
    switch print_type
     case 'ps'
      print('-dps',['plots/',name,'.ps'])
     case 'png'
      print('-dpng',['plots/',name,'.png']);
     case 'fig'
      saveas(fig_hl,['plots/',name,'.fig'],'fig');
    end
  end
end


plot_num = plot_num + 1;
if doplots(plot_num)
  fig_hl=figure(figno); figno = figno + 1;
  set(fig_hl,'Position', screen_position);
  screen_position = screen_position + offset;
  set(fig_hl,'PaperPosition',paper_position)

  subplot(3,1,1)
  plot(t,f(:,1))
  ylabel('f^T_1 (N)')
  axis([min(t) max(t) min(f(:,1))*1.1 max(f(:,1))*1.1])
  grid
  subplot(3,1,2)
  plot(t,f(:,2))
  ylabel('f^N_1 (N)')
  axis([min(t) max(t) min(f(:,2))*.9 max(f(:,2))*1.1])
  grid
  subplot(3,1,3)
  plot(t,f(:,1)./f(:,2))
  ylabel('f^T_1/f^N_1');
  xlabel('time (sec)')
  grid
  if doprints
    name = 'stance_leg_forces';
    switch print_type
     case 'ps'
      print('-dps',['plots/',name,'.ps'])
     case 'png'
      print('-dpng',['plots/',name,'.png']);
     case 'fig'
      saveas(fig_hl,['plots/',name,'.fig'],'fig');
    end
  end
end