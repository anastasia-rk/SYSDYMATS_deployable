local_init;
%% Test Van-der-Pol oscilator
% create folder to store simulation data
folderName = make_folder('../SYSDYMATS_data/vdpo');                     
% Initial conditions
tspan = [0 1000];                                                           % integration span
y0 = [2; 0];                                                                % ODE initial conditions
phase(1) = 0;                                                               % initial phase
w(1) = 2*rand;                                                              % initial frequency
d = 50; ak = 0.2;                                                           % excitation signal parameters
extract = [1:300];                                                         % extract from the array for the plot
% mu_array = [0.0625 0.125 0.25 0.3 0.5 0.8 1];
mu_array = [0.1:0.1:1];
sampFr = 25; % 2 % 5 % 10 % 20 % 50 % 100 % 120 % 150 % 200
t0 = [tspan(1):1/sampFr:tspan(end)]';
nFr = 1000;
f = sampFr*(0:nFr-1)/nFr;
for iMu = 1:length(mu_array)
    for k=2:d
        w(k) = 2*rand;
        phase(k) = phase(1) - k*pi*(k-1)/d;
    end
%% Generate outputs
    Mu = mu_array(iMu);
    ode = @(t,y) vdpo(t,y,d,w,phase,ak,Mu);
    [t,y] = ode45(ode, tspan, y0);
    clear u
%% Display the excitation
    for it = 1:length(t)
        for k=1:d
            fun(k) = ak*cos(2*pi*w(k)*t(it) + phase(k)); 
        end
        u(it,1) = sum(fun);
    end
    
%% Save data to matfile
yy0 = interp1(t,y,t0);
u0 = interp1(t,u,t0);
U_jw = fft(u0,nFr);
Y_jw = fft(yy0,nFr);
m_u = abs(U_jw);                                % Magnitude
p_u = angle(U_jw);                              % Phase
m_y = abs(Y_jw);                                % Magnitude
p_y = angle(Y_jw);                              % Phase

fileData = [t0, u0, yy0];
fileName = [folderName,'/',num2str(iMu),'V',];
save(fileName, 'fileData');
%% Plot solution
figName = ['mu = ',num2str(Mu)];
fig(figName,visFlag);
subplot(3,2,1)
plot(t,u); hold on;
plot(t0(extract),u0(extract),'o'); hold on;
xlim([t0(extract(1)) t0(extract(end))]);
xlabel('$t$')
ylabel('$u(t)$')
subplot(3,2,2)
plot(t,y(:,1)); hold on;
plot(t,y(:,2)); hold on;
plot(t0(extract),yy0(extract,1),'o'); hold on;
plot(t0(extract),yy0(extract,2),'*');
xlim([t0(extract(1)) t0(extract(end))]);
legend('$y(t)$','$\dot{y}(t)$')
xlabel('$t$')
ylabel('$y(t)$')
subplot(3,2,3)
plot(f,m_u)
xlabel('f,Hz'); ylabel('$|U(j\omega)|$');
subplot(3,2,4)
plot(f,m_y(:,1))
xlabel('f,Hz'); ylabel('$|Y(j\omega)|$');
subplot(3,2,5)
plot(f,rad2deg(p_u))
xlabel('f,Hz'); ylabel('$\angle U(j\omega)$');
subplot(3,2,6)
plot(f,rad2deg(p_y(:,1)))
xlabel('f,Hz'); ylabel('$\angle Y(j\omega)$');
end

%% Save external parameters 
params = 'mu_vanderpol';
values = mu_array';
fileName = 'External_parameters_V';
save(fileName,'params','values','sampFr');