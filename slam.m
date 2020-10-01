%% I. Initialization
% motion model noise
q = [0.01;0.01];
Q = diag(q.^2);
% measurement noise
m = [.1; 0.1];
M = diag(m.^2);
% R: robot initial pose
% u: control
% *****************************************%
R = [0;-1;0];
u = [0.1;6/180*pi];
% *****************************************%
% set landmarks
Marks = landmarks();
% sensor radius
% *****************************************%
sensor_r = 1.5;
% *****************************************%
% i-th landmark is sensored before, LANDMARK(i) = 1;
% j-th landmark is not sensored, LANDMARK(j) = 0;
LANDMARK = zeros(1,size(Marks,2));
% y_news -- landmarks first seen 
% y_olds -- landmarks has been sensored before
y_olds = zeros(3,size(Marks,2));
y_news = zeros(3,size(Marks,2));
%   State and covariance intialization
y = zeros(numel(R)+numel(Marks), 1); %** State && Map **%
P = zeros(numel(y),numel(y));
signature = zeros(1,size(Marks,2));
r = [1 2 3];
y(r) = R;
sigma = 0;
% Map starts at 4&5
s = [4 5];
% 60 s/circle
 loop =120;
 % poses_ -- store the actual state
poses_ = zeros(3,loop);
% poses -- store the estimated state
poses = zeros(3,loop);
 %  PLOT setup
mapFig = figure(1);
cla;
axis([-2.5 2.5 -2.5 2.5])
axis square
% landmarks 
WG = line('parent',gca,...
    'linestyle','none',...
    'marker','*',...
    'color','m',...
    'xdata',Marks(1,:),...
    'ydata',Marks(2,:));
% actual path
RG = line('parent',gca,...
    'marker','+',...
    'color','r',...
    'xdata',R(1),...
    'ydata',R(2));
% estimated path
rG = line('parent',gca,...
    'linestyle','none',...
    'marker','+',...
    'color','b',...
    'xdata',y(r(1)),...
    'ydata',y(r(2)));
% estimated landmarks
lG = line('parent',gca,...
    'linestyle','none',...
    'marker','+',...
    'color','k',...
    'xdata',[],...
    'ydata',[]);
% estimated landmarks covariance
eG1 = zeros(1,size(Marks,2));
for i = 1:numel(eG1)
    eG1(i) = line(...
        'parent', gca,...
        'color','k',...
        'xdata',[],...
        'ydata',[]);
end
% robot covariance
reG = line(...
    'parent', gca,...
    'color','r',...
    'xdata',[],...
    'ydata',[]);
% sensor (the center is the actual position of robot)
sensor1 = line(...
    'parent', gca,...
    'color','m',...
    'xdata',[],...
    'ydata',[],...
    'LineStyle','--');
sensor2 = line(...
    'parent', gca,...
    'color','m',...
    'xdata',[],...
    'ydata',[],...
    'LineStyle','--');
 true_pose = line(...
    'parent', gca,...
    'color','r',...
    'xdata',[],...
    'ydata',[],...
    'LineWidth',0.8);
  estimate_pose = line(...
    'parent', gca,...
    'color','b',...
    'xdata',[],...
    'ydata',[],...
    'LineWidth',0.8);
  %% II. EKF-SLAM
for t = 1:loop
    % control noise
    n = q.*randn(2,1);
    % return the robot actual position£»
    R = move(R, u, n);
    
    i_olds=1;
    i_news=1;
    
    for i = 1:size(Marks,2)
        %mesurement error
        v = m.*randn(2,1);
         yi= project(R, Marks(:,i)) + v;
        if yi(1) < sensor_r && LANDMARK(i) == 1
               y_olds(:,i_olds) = [yi(1);yi(2);i];
               i_olds = i_olds + 1;
        elseif  yi(1) < sensor_r &&  LANDMARK(i) == 0 
                y_news(:,i_news) = [yi(1);yi(2);i];
                i_news = i_news + 1;
                LANDMARK(i) = 1;                
        end
    end
   
    for i = i_olds:size(Marks,2)
        y_olds(:,i) = [100;0;0];
    end
    for i = i_news:size(Marks,2)
        y_news(:,i) = [101;0;0];
    end
  
    % EKF
    %  prediction
    [y(r), JacobianMotion, Vt] = move(y(r), u, [0 0]);
    P_rr = sigma;
    P(r,:) = JacobianMotion*P(r,:);
    P(:,r) = P(r,:)';
    sigma = JacobianMotion*P_rr*JacobianMotion' + Vt*Q*Vt';
    
     %  update
    end_old = find(y_olds(1,:)==100,1);  
    if isempty(end_old)
        end_old=size(y_olds,2)+1;
    end
    
    for j = 1:(end_old-1)
        % expectation
        if isempty(j)
            break
        end
        id = find(signature==y_olds(3,j),1);
        v = [id*2+2 id*2+3];
        [e, E_r, E_l] = project(y(r), y(v));
        H = [E_r E_l];
        rl   = [r v];
        E    = H * P(rl,rl) * H';
        
        % measurement
        yi_1 = y_olds(:,j);
        yi1 = yi_1(1:2,1);
        
        % innovation
        z = yi1 - e;
        if z(2) > pi
            z(2) = z(2) - 2*pi;
        end
        if z(2) < -pi
            z(2) = z(2) + 2*pi;
        end
        S = E+M;
        
        % Kalman gain
        K = P(:, rl) * H' * S^-1;
        
        % update
        y = y + K * z;
        P = P - K * S * K';

    end
    
    % for the landmarks which are never seen before   
    end_new = find(y_news(1,:)==101,1);
    if isempty(end_new)
        end_new=size(y_news,2)+1;
    end
    for m1 = 1:(end_new-1)
        if isempty(m1)
            break
        end
        id = find(signature==0,1);
        signature(id) = y_news(3,m1);
        
        % measurement
        yi_2 = y_news(:,m1);
        yi2 = yi_2(1:2,1);
        [y(s), L_r, L_y] = backProject(y(r ), yi2);
        P(s,:) = L_r * P(r,:);
        P(:,s) = P(s,:)';
        P(s,s) = L_r * sigma * L_r' + L_y * M * L_y';
        s = s + [2 2];
    end
    
   
    % obtain states info
    % estimated
    poses(1,t) = y(1);
    poses(2,t) = y(2);
    poses(3,t) = y(3);   
    % actual
    poses_(1,t) = R(1);
    poses_(2,t) = R(2);
    poses_(3,t) = R(3);
    
     % 5. PLOT
     % Actual positon and the range of sensor
    set(RG, 'xdata', R(1), 'ydata', R(2));
    circle_x = linspace((R(1)-0.9999*sensor_r),(R(1)+0.9999*sensor_r));
    circle_y1 = sqrt(sensor_r^2 - (circle_x - R(1)).^2) + R(2);
    circle_y2 = R(2) - sqrt(sensor_r^2 - (circle_x - R(1)).^2);
    set(sensor1,'xdata',circle_x,'ydata',circle_y1);
    set(sensor2,'xdata',circle_x,'ydata',circle_y2);
    
    % Estimated position
    set(rG, 'xdata', y(r(1)), 'ydata', y(r(2)));
    Circle_x = linspace((y(r(1))-0.9999*sensor_r),(y(r(1))+0.9999*sensor_r));
    Circle_y1 = sqrt(sensor_r^2 - (Circle_x - y(r(1))).^2) + y(r(2));
    Circle_y2 = y(r(2)) - sqrt(sensor_r^2 - (Circle_x - y(r(1))).^2);

    
    % Robot actual&estimated path
    set(estimate_pose,'xdata',poses(1,1:t),'ydata',poses(2,1:t));
    set(true_pose,'xdata',poses_(1,1:t),'ydata',poses_(2,1:t));
    
    legend([estimate_pose true_pose lG WG],{'Estimated Path','Actual Path' 'Estimated landmark' 'Actual landmark'})

  
  if s(1)==4
        continue
  end
 
  % The estimated locations of landmark
  w = 2:((s(1)-2)/2);
  w = 2*w;
  landmark_estimated_x = y(w);
  landmark_estimated_y = y(w+1);
  set(lG, 'xdata', landmark_estimated_x, 'ydata', landmark_estimated_y);
  
%%%%% 1- the landmark is never sensored before (BLUE)
  for g1 = 1:(end_new-1)
      if isempty(g1)
            break
      end
      o1 = y_news(3,g1);
      h1 = find(signature==o1,1);
      temp1 = [2*h1+2;2*h1+3];
      le = y(temp1);
      LE = P(temp1,temp1);
      [X,Y] = cov2elli(le,LE,3,16);   
      set(eG1(o1),'xdata',X,'ydata',Y,'color','b');
  end
  %%%% 2- the landmark has been sensored and is sensored again (RED)
  for g2 = 1:(end_old-1)
      if isempty(g2)
            break
      end
      o2 = y_olds(3,g2);
      h2 = find(signature==o2,1);
      temp2 = [2*h2+2;2*h2+3];
      le = y(temp2);
      LE = P(temp2,temp2);
      [X,Y] = cov2elli(le,LE,3,16);  
      set(eG1(o2),'xdata',X,'ydata',Y,'color','r');
  end
  %%%% 3- the landmark has been sensored and is NOT sensored now (BLACK)
  v = find(signature==0,1);
  if isempty(v)
      v = size(signature,2)+1;
  end
  for g3 = 1:v-1
      if isempty(g3)
            break
      end
      a = find(y_olds(3,:)==signature(g3),1);
      b = find(y_news(3,:)==signature(g3),1);   
      if (isempty (a)) && (isempty(b)) 
         temp3 =  [2*g3+2;2*g3+3];
            le = y(temp3);
      LE = P(temp3,temp3);
      [X,Y] = cov2elli(le,LE,3,16);
      set(eG1(signature(g3)),'xdata',X,'ydata',Y,'color','k');
      end
  end

% Covariance of the robot position (RED)
%      if t > 1
%          re = y(r(1:2));
%          RE = P(r(1:2),r(1:2));
%          [X,Y] = cov2elli(re,RE,3,16);
%          set(reG,'xdata',X,'ydata',Y);
%      end
%    
   drawnow;
   
   pause(0.5);
end   