%% Figure5

% Author: Jesse Sharp; Last Update: 21/09/2018

%

% Corresponds to Figure 5



%% Set-up
clear

Tfinal = 50;  %Specified final time
dt = 0.001; %Time-step
N = Tfinal/dt+1; %Number of nodes in time discretisation
t_y = linspace(0,Tfinal,N); %Time discretisation for plots
omega = 0.9; %Portion of previous iteration's control maintained when updating control
RelTol = 1e-3; %Desired relative tolerance for convergence
MaxIters = 250; %Number of iterations to perform before giving up if convergence is not reached
U = zeros(1,N); %Initial guess for the control

%Model parameters
ps = 0.5; %Proliferation of S
pa = 0.43; %Proliferation of A
pl = 0.27; %Proliferation of L
gs = 0.14; %Differentiation of S to A
ga = 0.44; %Differentiation of A to D
gl = 0.05; %Differentiation of L to T
ux = 0.275; %Migration of D into the blood stream
ut = 0.3; %Migration of T into the blood stream
K1 = 1; %Carrying capacity of the compartment with S
K2 = 1; %Carrying capacity of the compartment with SA + L
Alpha = 0.015; %Michaelis-Menten kinetic parameter alpha
Gamma = 0.1; %Michaelis-Menten kinetic parameter gamma

%Pay-off weightings
a1 = 1; %Weighting on negative impact of control
a2 = 1; %Weighting on negative impact of leukaemia

%initialisation using the approximate steady state of the state system
%for the given parameters
S(1) = 1-gs/ps; 
A(1) = 0.3255; 
X(1) = 0.5207; 
L(1) = 0.3715; 
T(1) = 0.0619; 
y(:,1) = [S(1),A(1),X(1),L(1),T(1)];

%Other set-up
iterations = 0; %initialise iteration count
RelTolStore = []; %Initialise vector to monitor convergence

%State equations
State = @(t,y,U) [ps*y(1)*(K1-y(1))-gs*y(1);
    gs*y(1)+pa*y(2)*(K2-(y(2)+y(4)))-ga*y(2);
    ga*y(2)-ux*y(3);
    pl*y(4)*(K2-(y(2)+y(4)))-gl*y(4)-Alpha*y(4)/(Gamma+y(4))-U*y(4);
    gl*y(4)-ut*y(5)];

%Costate equations
Costate = @(t,Lambda,y,U) [-(-2*y(1)*Lambda(1)*ps-gs*Lambda(1)+gs*Lambda(2)+Lambda(1)*ps);
    -(-2*y(2)*Lambda(2)*pa-y(4)*Lambda(2)*pa-y(4)*Lambda(4)*pl-ga*Lambda(2)+Lambda(2)*pa+ga*Lambda(3));
    -(-ux*Lambda(3));
    -(2*a2*y(4)-pa*y(2)*Lambda(2)-Lambda(4)*pl*y(2)-2*pl*y(4)*Lambda(4)+Lambda(4)*pl-Lambda(4)*gl-Lambda(4)*Alpha/(Gamma+y(4))+Lambda(4)*Alpha*y(4)/(Gamma+y(4))^2-U*Lambda(4)+gl*Lambda(5));
    ut*Lambda(5)];

 %% Forward-backward sweep
while(iterations<MaxIters)
    
uold=U; %Store control from previous iteration

i = 0; %Initialise loop variable
t = 0; %Initialise time for forward sweep

%Forward sweep using fourth-order Runge-Kutta scheme
while i < N-1
    t = t + dt;
    i=i+1;
    k1 = State(t,y(:,i),U(i));
    k2 = State(t,y(:,i)+dt*k1/2,0.5*(U(i)+U(i+1)));
    k3 = State(t,y(:,i)+dt*k2/2,0.5*(U(i)+U(i+1)));
    k4 = State(t,y(:,i)+dt*k3,U(i+1));
    y(:,i+1) = y(:,i) + (dt/6)*(k1+2*k2+2*k3+k4);
end

t = Tfinal; %Initialise time for backward sweep
Lambda = zeros(5,length(y)); %Initialise Lambda
Lambda(:,end) = [0,0,0,0,0]; %Apply transversality conditions to obtain final time condition on Lambda (costate)
i = 0; %Initialise loop variable
j = N; %Initialise loop variable

%Backward sweep using fourth-order Runge-Kutta scheme 
while j > 1 
    t = t - dt;
    i = i+1;
    j = N-i;
    k1 = Costate(t,Lambda(:,j+1),y(:,j+1),U(j+1));
    k2 = Costate(t,Lambda(:,j+1)-dt*k1/2,0.5*(y(:,j)+y(:,j+1)),0.5*(U(j)+U(j+1)));
    k3 = Costate(t,Lambda(:,j+1)-dt*k2/2,0.5*(y(:,j)+y(:,j+1)),0.5*(U(j)+U(j+1)));
    k4 = Costate(t,Lambda(:,j+1)-dt*k3,y(:,j),U(j));
    Lambda(:,j) = Lambda(:,j+1) - (dt/6)*(k1+2*k2+2*k3+k4);
end

%Update the value of the control and apply a relaxation factor.
Uupdate = y(4,:).*Lambda(4,:)/(2*a1); %Calculated updated control control
U = omega*U + (1-omega)*Uupdate; %Actual updated control after applying relaxation to aid convergence

% %To view the solution as it converges, uncomment this block (line 111 to line 126) to produce interim
% figures
% box on
% line1 = plot(t_y,y(1,:),'LineWidth',2);
% hold on
% line2 = plot(t_y,y(2,:),'LineWidth',2);
% line3 = plot(t_y,y(3,:),'LineWidth',2);
% line4 = plot(t_y,y(4,:),'LineWidth',2);
% line5 = plot(t_y,y(5,:),'LineWidth',2);
% line6 = plot(t_y,U,'k--','LineWidth',2);
% hL = legend([line1,line2,line3,line4,line5,line6],{'S','A','D','L','T','u*'},'Location','northeast');
% ylabel('State','fontsize',18);
% xlabel('Time','fontsize',18);
% axis([0,Tfinal,0,1])
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 18)
% hold off
% drawnow

%Check for convergence
RelTolTest = RelTol*norm(U) - norm(U-uold);
RelTolStore = [RelTolStore RelTolTest];
if RelTolTest > 0 
    fprintf('Specified relative tolerance of %g has been met \n\r',RelTol)
    break
end  

%Update and display iteration count
iterations = iterations+1;
fprintf('Iterations completed: %d  \n\r',iterations) 
end


%% Plot final state with optimal control
colours = [ 
    222/255  125/255  0/255 
    237/255  177/255  32/255 
    255/255  255/255  0/255
    20/255  42/255  140/255
    0/255  114/255  189/255
]; %Define colours for plot
 
figure
set(gca, 'ColorOrder', colours);
hold on
box on
line1 = plot(t_y,y(1,:),'LineWidth',2);
line2 = plot(t_y,y(2,:),'LineWidth',2);
line3 = plot(t_y,y(3,:),'LineWidth',2);
line4 = plot(t_y,y(4,:),'LineWidth',2);
line5 = plot(t_y,y(5,:),'LineWidth',2);
line6 = plot(t_y,U,'k--','LineWidth',2);
hL = legend([line1,line2,line3,line4,line5,line6],{'S','A','D','L','T','u*'},'Location','northeast');
ylabel('State','fontsize',18);
xlabel('Time','fontsize',18);
axis([0,Tfinal,0,1])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 18)
