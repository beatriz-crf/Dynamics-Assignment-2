%% Load Structure 

close all; clc; clear;

%Loadstructure 
[~,xy,nnod,~,idb,ndof,incid,l,gamma,m,EA,EJ,posit,~,pr]=loadstructure ;

% draw the structure
dis_stru(posit,l,gamma,xy,pr,idb,ndof);

% assemble mass and stiffness matrices
[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

%% Verify Lmax
omega0_element = zeros(length(l),1);  % first nat frequency of each element
fmax = 20; %Hz
omax = fmax*2*pi; %rad/s
SF = 2; % safety factor

% Calculate permissible Lmax of each beam type
EJ_red = 5.7380e6;  % [Nm^2]
m_red = 26.2;   % mass per unit length [kg/m]
Lmax_red = sqrt(pi^2/(SF*omax)*sqrt(EJ_red/m_red));

EJ_green = 3.5400e5;  % [Nm^2]
m_green = 8.1;   % mass per unit length [kg/m]
Lmax_green = sqrt(pi^2/(SF*omax)*sqrt(EJ_green/m_green));

EJ_blue = 1.7995e6;  % [Nm^2]
m_blue = 15.8;   % mass per unit length [kg/m]
Lmax_blue = sqrt(pi^2/(SF*omax)*sqrt(EJ_blue/m_blue));

% Verify appropriateness of mesh by element length criterion
for i=1:length(l)
    omega0_element(i) = ((pi/l(i))^2)*sqrt((EJ(i)/(m(i))));
    
    if omega0_element(i)/omax <= (SF) % safety condition to verify Lmax
 fprintf ('safety factor is less than 2.0 and therefore not ok \n');
    end
end

%% Mode shape calculation

% Eigenvector problem
MFF = M(1:ndof, 1:ndof);    % free-free mass matrix
KFF = K(1:ndof, 1:ndof);    % free-free stiffness matrix

[modes, omega0_square] = eig(MFF\KFF);

omega0 = sqrt(diag(omega0_square));

[omega0,i_omega0] = sort(omega0);  % Sort frequencies in ascending order 
freq0 = omega0/2/pi;
modes = modes(:,i_omega0);   % Sort mode shapes in ascending order 

% Mode shape plot
scale_factor = 3;
nmodes = 3;

for i = 1:nmodes
    figure;
    diseg2(modes(:,i),scale_factor,incid,l,gamma,posit,idb,xy);
    title(sprintf('Mode %d -- f_o_%d = %.3f [Hz]',i,i,freq0(i)))
end

%% Frequency Response Functions
alpha = 0.1;
beta = 2.0e-4;

C = alpha*(M) + beta*(K);   % damping matrix
CFF = C(1:ndof, 1:ndof);    % free-free damping matrix

freq_res = 0.01;    % frequency resolution [Hz]
freq = 0:freq_res:fmax;   % frequency vector [Hz]
Omega = 2*pi*freq;  % frequency vector [rad/s]

% Build force vector F0
% Vertical input force applied at point A (node 3)
F0 = zeros(ndof,1);
F0_idx = idb(3,2);
F0(F0_idx) = 1;

% Build output response X0
% Vertical displacement at point A (node 3)
% Vertical displacement at point B (node 11)
X0 = zeros(ndof,length(freq));
X0_idx1 = idb(3,2);     % output measurement location
X0_idx2 = idb(11,2);    % output measurement location

for ii = 1:length(Omega)
    A = -Omega(ii)^2*MFF + 1i*Omega(ii)*CFF + KFF;  % inverse of FRF
    X0(:,ii) = A\F0;    % output response
end

% FRF: G=output/input
G1 = X0(X0_idx1,:)/F0(F0_idx); % FRF in correspondence of output location 1
G2 = X0(X0_idx2,:)/F0(F0_idx); % FRF in correspondence of output location 2

% Plot FRF
figure; 
subplot(2,1,1); 
semilogy(freq, abs(G1),'LineWidth',2,'DisplayName','Point A')

set(groot, 'defaultTextInterpreter', 'tex'); 
title('FRF: Outputs at Points A & B; Input at Point A');
xlabel('Frequency [Hz]'); ylabel('|G|')
legend('show', 'Location', 'northeast'); grid on;
% axis([0 200 10e-6 5e-1]);

% Plot Exp FRF Phase
subplot(2,1,2);
plot(freq, angle(G1), 'LineWidth', 2, 'DisplayName', 'Point A');
grid minor; 
xlabel('Frequency [Hz]');
ylabel('∠G [rad]');
legend('show', 'Location', 'northeast');

%% Modal Analysis


%%
