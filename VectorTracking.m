clc
clear all
Hs=60;%m
select = 'B';
a=0;b=0;c=0;d=0;p=0;
switch(select)
    case 'B'
        a = 0.310;
        b = 0.897;
        c = 0.1393;
        d = 0.9467;
        p =0.15;
    case 'D'          %DD = 4; 2<u<=5
        a = 0.122;
        b = 0.916;
        c = 0.0856;
        d = 0.8650;
        p =0.3;
    case 'E'          %E = 6;   u<=2
        a = 0.0934;
        b = 0.912;
        c = 0.1094;
        d = 0.7657;
        p =0.3;
end

z0 = 10;
u0=2.0;%m/s wind velocity
Qs = [4419.834747
4419.810703
3665.279742
3680.133091
3705.051432
3730.263523
3750.461482]; %g



vg = 1.0; %m/s gas exit velocity
dm = 2.0; %m stack diameter
t0 = 333.15; %K gas exit temperature 60C
t1 = 293.15; %K ambient temperature
f0 = 3.12*0.785*vg*dm^2*(t0-t1)/t0;
if (f0 > 55)
    x0 = 34*exp(0.4*log(f0));
else
    x0 = 14*exp(0.625*log(f0));
end
He = Hs + 1.6*exp(log(f0)/3)*exp(2*log(3.5*x0)/3)/u0;
X=-20:1:100;
Y=-20:1:20;
Z=[74];

dt=1;

vx=-6; %m/s
vy=0; %m/s
vessel_position = [0,0];
emission_positions = [vessel_position;
    vessel_position+[-vx*dt,-vy*dt];
    vessel_position+[-vx*2*dt,-vy*2*dt];
    vessel_position+[-vx*3*dt,-vy*3*dt];
    vessel_position+[-vx*4*dt,-vy*4*dt];
    vessel_position+[-vx*5*dt,-vy*5*dt];
    vessel_position+[-vx*6*dt,-vy*6*dt]];
C_max = zeros(size(Z));
C_area = zeros(size(Z));
for h=1:size(Z,2)
    z = Z(h);
    uz = u0*(z/z0)^p;
    emissions_num = size(emission_positions,1);
    sum_C = zeros(size(X,2),size(Y,2));
    for n=1:emissions_num
        
        C = zeros(size(X,2),size(Y,2));
        for i = 1:size(X,2)
            for j = 1:size(Y,2)
        
                x = X(i)-emission_positions(n,1);
                y = Y(j)-emission_positions(n,2);
                if(x>0)
                    xx=abs(x);
                    yy=abs(y);
                    sigma_x=a*(xx.^b);
                    sigma_y=sigma_x;%sigma_x;
                    sigma_z=c*xx.^d;
                    K_0=(2*pi).^(1.5)*sigma_y*sigma_z*sigma_y;
                    K_1 = 1/K_0;        
                    K_2= -(xx-uz*n*dt).^2/(2*sigma_x.^2);
                    K_3= -(yy).^2/(2*sigma_y.^2);
                    K_4= -(z-He).^2/(2*sigma_z.^2);
                    K_5= -(z+He).^2/(2*sigma_z.^2);
                    C(i,j)=(Qs(n)*K_1*exp(K_2)*exp(K_3)*(exp(K_4)+exp(K_5)));
                    %C(i,j)=integral(@(tba)(Qs(n)*K_1*exp((xx-uz*tba).^2*K_2)*exp(K_3)*(exp(K_4)+exp(K_5))),0,dt*n);
                    C(i,j)=real(C(i,j));
                else
                    C(i,j) = 0;
                end
        
            end
        end
        result_C = findNaN(C);
        sum_C = sum_C+result_C;
    end
    
    for num = 1:20
        C_noise = zeros(size(X,2),size(Y,2));
        mu = [round(rand*60), round(rand*20-10)];
        Sigma_1 =3;
        Sigma_2 =3;
        for i = 1:size(X,2)
            for j = 1:size(Y,2)
                C_noise(i,j) = exp(-0.5*(mu(1)-X(i)).^2/Sigma_1.^2-0.5*(mu(2)-Y(j)).^2/Sigma_2.^2)/(2*pi*Sigma_1*Sigma_2)*sum_C(20+mu(1),20+mu(2));
            end
        end
        sum_C = sum_C+C_noise;
        
    end    
    C_max(h) = max(max(sum_C));
    [rows,cols] = find(sum_C==C_max);
    
end
n = size(rows,1);
ais_pos = [sum(rows),sum(cols)]/n+[-20,-20];
ais_mea = C_max(1);
init_pos = [80,-5];
pos = init_pos;
pos_opt = init_pos;
pos_opt_v = init_pos;
velocity = 1; %m/s
dire = [-1,0];
dy_velocity = 1;

pos_record(1,:) = init_pos;
pos_opt_record(1,:) = pos_opt;
pos_opt_v_record(1,:) = pos_opt_v;

%%only tracking
mea_buf = [];pos_buf = [];t = 1;
while(t<=200)
    mea = sum_C(20+round(pos(1)),20+round(pos(2)));
    
    mea_record(t,:) = mea;
    
    if mea >= ais_mea*0.9
        fprintf('only tracking end: %f,%f\n',pos(1),pos(2))
        break;
    end
    if size(mea_buf,1)>=5
       mea_buf(1)=[];
       pos_buf(1,:)=[];
       mea_buf(10) = mea;
       pos_buf(10,:) = pos;
    else
        mea_buf(size(mea_buf,1)+1,:) = mea;
        pos_buf(size(pos_buf,1)+1,:) = pos;
        disp(mea_buf);
    end
    vector(t,:) = [0,0];
    for i=1:size(mea_buf,1)-1
        vector1 = pos_buf(size(pos_buf,1),:)-pos_buf(i,:);
%         fprintf('vector1: %f,%f\n',vector1(1),vector1(2))
        
        if abs(vector1(1))>0 || abs(vector1(2))>0
            vector1 = vector1/norm(vector1,2);
            diff_c = mea_buf(size(mea_buf,1),:)-mea_buf(i,:);
            vector1 = diff_c*vector1;
%             fprintf('vector1: %f,%f\n',vector1(1),vector1(2))
            
        end
        vector(t,:) = vector(t,:) +vector1;
    end
    
    if vector(t,1)>0 || vector(t,2)>0
        vector(t,:) = vector(t,:)/norm(vector(t,:),2);        
    else
        vector(t,:) = dire;
    end
    pos = pos + velocity*vector(t,:);
    t= t+1;
    pos_record(t,:) = pos;
    if pos(1)<-20 || pos(2)<-20 || pos(1)>=100 || pos(2)>=20
        fprintf('failed: %f,%f\n',pos(1),pos(2))
        break;
    end
end
% optimal
t = 1;mea_buf = [];pos_buf = [];
while(t<=200)
  
    mea_opt = sum_C(20+round(pos_opt(1)),20+round(pos_opt(2)));
    mea_opt_record(t,:) = mea_opt;
    if mea_opt >= ais_mea*0.9
        fprintf('opt end: %f,%f time cost: %f\n',pos_opt(1),pos_opt(2),t)
        break;
    end
    if size(mea_buf,1)>=5
       mea_buf(1)=[];
       pos_buf(1,:)=[];
       mea_buf(10) = mea_opt;
       pos_buf(10,:) = pos_opt;
    else
        mea_buf(size(mea_buf,1)+1,:) = mea_opt;
        pos_buf(size(pos_buf,1)+1,:) = pos_opt;
        disp(mea_buf);
    end
    vector(t,:) = [0,0];
    for i=1:size(mea_buf,1)-1
        vector1 = pos_buf(size(pos_buf,1),:)-pos_buf(i,:);
%         fprintf('vector1: %f,%f\n',vector1(1),vector1(2))
        
        if abs(vector1(1))>0 || abs(vector1(2))>0
            vector1 = vector1/norm(vector1,2);
            diff_c = mea_buf(size(mea_buf,1),:)-mea_buf(i,:);
            vector1 = diff_c*vector1;
%             fprintf('vector1: %f,%f\n',vector1(1),vector1(2))
            
        end
        vector(t,:) = vector(t,:) +vector1;
    end
    
    if vector(t,1)>0 || vector(t,2)>0
        vector(t,:) = vector(t,:)/norm(vector(t,:),2);        
    else
        vector(t,:) = dire;
    end
 
    vector2 = ais_pos-pos_opt;
    vector2 = vector2/norm(vector2,2);
    diff_c = ais_mea-mea_opt;
    vector2 = diff_c*vector2;
    vector_opt(t,:) = vector(t,:)+vector2;
    if vector_opt(t,1)>0 || vector_opt(t,2)>0
        vector_opt(t,:) = vector_opt(t,:)/norm(vector_opt(t,:),2);
    else
        vector_opt(t,:) = dire;
        
    end

    pos_opt = pos_opt + velocity*vector_opt(t,:);
    t= t+1;
    pos_opt_record(t,:) = pos_opt;

    if pos_opt(1)<-20 || pos_opt(2)<-20 || pos_opt(1)>=100 || pos_opt(2)>=20
        fprintf('opt failed: %f,%f\n',pos_opt(1),pos_opt(2))
        break;
    end
end


%%%%update
t = 1;mea_buf = [];pos_buf = [];
while(t<=200)
    mea_opt_v = sum_C(20+round(pos_opt_v(1)),20+round(pos_opt_v(2)));
    mea_opt_v_record(t,:) = mea_opt_v;
    if mea_opt_v >= ais_mea*0.9
        fprintf('opt_v end: %f,%f; time cost: %f',pos_opt_v(1),pos_opt_v(2),t)
        break;
    end
    if size(mea_buf,1)>=5
       mea_buf(1)=[];
       pos_buf(1,:)=[];
       mea_buf(10) = mea_opt_v;
       pos_buf(10,:) = pos_opt_v;
    else
        mea_buf(size(mea_buf,1)+1,:) = mea_opt_v;
        pos_buf(size(pos_buf,1)+1,:) = pos_opt_v;
        disp(mea_buf);
    end
    vector(t,:) = [0,0];
    for i=1:size(mea_buf,1)-1
        vector1 = pos_buf(size(pos_buf,1),:)-pos_buf(i,:);
%         fprintf('vector1: %f,%f\n',vector1(1),vector1(2))
        
        if abs(vector1(1))>0 || abs(vector1(2))>0
            vector1 = vector1/norm(vector1,2);
            diff_c = mea_buf(size(mea_buf,1),:)-mea_buf(i,:);
            vector1 = diff_c*vector1;
%             fprintf('vector1: %f,%f\n',vector1(1),vector1(2))
        end
        vector(t,:) = vector(t,:) +vector1;
    end
    if vector(t,1)>0 || vector(t,2)>0
        vector(t,:) = vector(t,:)/norm(vector(t,:),2);        
    else
        vector(t,:) = dire;
    end
 
    vector2 = ais_pos-pos_opt_v;
    vector2 = vector2/norm(vector2,2);
    diff_c = ais_mea-mea_opt_v;
    vector2 = diff_c*vector2;
    vector_opt(t,:) = vector(t,:)+vector2;
    if vector_opt(t,1)>0 || vector_opt(t,2)>0
        dy_velocity = norm(vector_opt(t,:),2);
        vector_opt(t,:) = vector_opt(t,:)/norm(vector_opt(t,:),2);
    
        if dy_velocity > 5
            dy_velocity = 5;
        end
%         fprintf('dy_velocity: %f\n',dy_velocity)
    else
        vector_opt(t,:) = dire;
        dy_velocity = 1;
    end
%     fprintf('pos_opt_v,dy_velocity: %f %f %f\n',pos_opt_v(1), pos_opt_v(2), dy_velocity)
    pos_opt_v = pos_opt_v + dy_velocity*vector_opt(t,:);
    t= t+1;
    pos_opt_v_record(t,:) = pos_opt_v;
%     fprintf('pos_opt_v_record: %f %f %f\n',pos_opt_v_record(t,1), pos_opt_v_record(t,2),t)
    if pos_opt_v(1)<=-20 || pos_opt_v(2)<=-20 || pos_opt_v(1)>=100 || pos_opt_v(2)>=20
        fprintf('opt_v failed: %f,%f\n',pos_opt_v(1),pos_opt_v(2))
        break;
    end
end


figure(2)
[F,h0] =contourf(X,Y,sum_C',50);
set(h0,'LineColor','none')
hold on
plot(emission_positions(1,1),emission_positions(1,2),'kp','markersize',15);
plot(pos_record(:,1),pos_record(:,2),'k.-','markersize',15);
plot(pos_opt_record(:,1),pos_opt_record(:,2),'y.-','markersize',15);
plot(pos_opt_v_record(:,1),pos_opt_v_record(:,2),'r.-','markersize',15);
%Make colourbar easier to read
h1=colorbar;
    %caxis([0,5])
colormap('jet');
set(h1,'fontsize',14);
ylabel(h1,'Concentration (g/m^3)');
%Change title based on stability condition provided
xlabel('X-position(m)','FontSize',14);
ylabel('Y-position(m)','FontSize',14);
hold off
figure(1)
subplot(2,1,1)
hold on
plot(pos_record(:,1),'k--');
plot(pos_record(:,2),'k--');
plot(pos_opt_record(:,1),'y-');
plot(pos_opt_record(:,2),'y-');
plot(pos_opt_v_record(:,1),'g-');
plot(pos_opt_v_record(:,2),'g-');
ylabel('position coordinates (m)','FontSize',14);
%Change title based on stability condition provided
xlabel('Time(s)','FontSize',14);
hold off
subplot(2,1,2)
hold on
plot(mea_record,'k--');
plot(mea_opt_record(:,1),'y-');
plot(mea_opt_v_record(:,1),'g-');
ylabel('gas measurements (g/m^3)','FontSize',14);
xlabel('Time(s)','FontSize',14);
hold off
