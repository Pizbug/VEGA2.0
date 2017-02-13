clc
clear all

%% PSI_RZ

fread=fopen('PSI_RZ.dat','r');
% fread = fopen('rzpsi_Circle.dat','r');

rz=fscanf(fread,'ZONE I=%d J=%d K=1\n');

rnum = rz(2);
znum = rz(1);

rzpsi = fscanf(fread,'%lf\t%lf\t%lf',[3,rnum*znum]);

%Z ¿˙¿Â
for j = 1 : znum

    z(j) = rzpsi(2,j);
    
end

i = 1;
j = 1;

for k = 1 : rnum*znum
 
    psi(i,j) = rzpsi(3,k);
    
    if( mod(k,znum) == 0 )
        
        r(i) = rzpsi(1,k);

        i = i + 1;
        j = 1;
        continue;
    
    end
    
    j = j + 1;
   
end

fclose(fread);

figure(11)
hold on
contour(r,z,psi', 100);


% plot(limx, limy)
% %%%%%%%%%%%%%
% %Plot option%
%%%%%%%%%%%%%
axis equal
% axis([1,2.4,-1.5,1.5])    %for KSTAR
% axis([0.07,0.9,-0.8,0.8]) % for VEST
xlabel('R','fontsize',15)
ylabel('Z','fontsize',15)
set(gca,'fontsize',15);
grid on
grid minor

% figure(2)
% rp = linspace(r(1,1), r(1,rnum), 105);
% zp = linspace(z(1,1), z(1,znum), 105);
% psip = interp2(r, z, psi, rp, zp, 'spline');
% 
% contour(rp,zp,psip, 50);
% % %%%%%%%%%%%%%
% % %Plot option%
% %%%%%%%%%%%%%
% axis equal
% % axis([1,2.4,-1.5,1.5])
% xlabel('R','fontsize',15)
% ylabel('Z','fontsize',15)
% set(gca,'fontsize',15);
% grid on
% grid minor
% 
% tilefigs

%%

BT = dlmread('BT.txt');

figure(22)
plot(BT(:,1), BT(:,2));
hold on





