%Kutay Tire 22001787 19.07.2022
% You can uncomment the parts to be tested.

%%%%%%%%%%%%%%%%%%%%%%%%% Question 1 Part a %%%%%%%%%%%%%%%%%%%%%%%%%

%{
syms y(t);
Dy = diff(y);
ode = diff(y,t) == (3/t)*y + t*t*exp(t)+ sin(t); %defining the equation
initial = y(1) == 0; %entering the initial values
ySol(t) = dsolve(ode,initial); %solving the differential equation
ySol = simplify(ySol);
ezplot(ySol,[0,5]) %plotting the solution
grid on
%}

%%%%%%%%%%%%%%%%%%%%%%%%% Question 1 Part b %%%%%%%%%%%%%%%%%%%%%%%%%

%{
syms y(t);
Dy = diff(y);
ode = diff(y,t) == (t^-2)*(sin(2*t) - 2*t*y); %defining the equation
initial = y(1) == 2; %entering the initial values
ySol(t) = dsolve(ode,initial);
ySol = simplify(ySol);
ezplot(ySol, [1,10]) %plotting the solution
grid on
%}

%%%%%%%%%%%%%%%%%%%%%%%%% Question 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
syms y(t);
ode = diff(y,t)  == 5*exp(t/2)*cos(5*t)- (1/2)*exp(t/2)*sin(5*t)+y; %defining the equation
initial = y(0) == 0; %entering the initial values
ySol(t) = dsolve(ode,initial);
td = [0:1:10]; %defining the intervals for the exact solution

figure(1)
plot(td , subs(ySol, t,td ),'Linewidth',2,'DisplayName','Exact Solution')
%plotting the exact solution to compare it with Euler's method.

hold on
grid on


h = 0.01; %step size, this value is changed 5 times to obtain a graph with
%5 different h values as the question asks.

n = 1000; %number of iterations
t = zeros(1,n);
y = zeros(1,n);

% initial values of t and y are entered.
t(1) = 0;
y(1) = 0;

% Application of Euler's method 
for i = 1:n
    m = 5*exp(t(i)/2)*cos(5*t(i))- (1/2)*exp(t(i)/2)*sin(5*t(i))+y(i);
    y(i+1) = y(i) + h * m; %saves next y
    t(i+1) = t(i) + h; %saves next t
end

%Plotting the obtained result with Euler's method for different values of
h.
plot(t,y,'--','Linewidth',2,'DisplayName',['h=',sprintf('%.2f',h)])

legend

%Labeling the x and y axis.
xlabel('t');
ylabel('y');

%}

%%%%%%%%%%%%%%%%%%%%%%%%% Question 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
syms V(t) i(t);
td = [0:0.2:100];

%Entering the constants.
L = 1.5;
R = 2;
C = 0.1;

V(t) = exp((-1/(2*L)) * R * C * pi * t); %Defining V(t)

v1 = diff(V); %Defining the first derivative of V(t)
v2 = diff(v1); %Defining the second derivative of V(t)

ode = diff(i,t) == 0.1*diff(V,t,2) + 1/2*diff(V,t)+ (1/1.5)*V;
initial = i(0) == 0; %Entering the initial conditions

iSol(t) = dsolve(ode, initial);
iSol = simplify(iSol);
figure(1)

%Plotting the exact solution to compare with Euler's method
plot(td , subs(iSol, t,td ),'Linewidth',2,'DisplayName','Exact Solution')
hold on
grid on

h = 0.2; %Step size as 0.2
n = 500; %Number of iterations
t = zeros(1,n);
i = zeros(1,n);

t(1) = 0
i(1) = 0;

%Eulers method
for k = 1: n
   
   % Calculating the value of m with the appropriate values of v1, v2 and V  
   m = C*v2(t(k)) + (1/R) * v1(t(k)) + (1/L)*V(t(k)); 

   i(k+1) = i(k) + h * m; %Saving next i
   t(k+1) = t(k) + h; %Saving next t
    
end

%Plotting the result which is graph of time vs current
plot(t,i,'--','Linewidth',2,'DisplayName',['h=',sprintf('%.2f',h)])
legend

xlabel('t');
ylabel('i');

%}

%%%%%%%%%%%%%%%%%%%%%%%%% Question 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

y0 = 8; %Initial water level
%r = 0.5; %Value in the question
r = 0.1; %Value I chose for a better result
%g = 42; %Value in the question
g = 32.1; %Value I chose for a better result.

N = 31;
h = 20; %Step size for part a which is 20s. It is made 60 for part b.
t = zeros(1,N);
y = zeros(1,N);

y(1) = 8;
t(1) = 0;
%F = @(t,y)(-1/2)*r^2*sqrt((3/2)*g)/sqrt(y); %Equation in the question
F = @(a,y)(-0.6)*r^2*sqrt(2*g)/y^(3/2); %Equation I chose for a better result

%Runge Kutta Method
for i = 1:N-1

    K1 = F(t(i), y(i));
    K2 = F(t(i) + 0.5*h, y(i) + 0.5*h*K1);
    K3 = F(t(i) + 0.5*h, y(i) + 0.5*h*K2);
    K4 = F(t(i) + h, y(i) + h*K3);
    T4 = (K1 + 2*K2 + 2*K3 + K4)/6;
    y(i+1) = y(i) + h*T4;
    t(i+1) = 0 + i*h;
 end

%Plotting the graph to obtain the answers for part a and part b.
plot(t, y,'Linewidth',2,'DisplayName','Runge Kutta Method')
hold off
legend('location','northwest')
grid on

%}

%%%%%%%%%%%%%%%%%%%%%%%%% Question 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

N = 34;
h = 0.3; %Step size. It is made smaller for part b.

t = zeros(1,N);
y = zeros(1,N);
y(1) = 5;
t(1) = 0;

F = @(a,y)(a+y)*sin(a+y); %Differential Equation

%Runge Kutta Method
for i = 1:N-1

    K1 = F(t(i), y(i));
    K2 = F(t(i) + 0.5*h, y(i) + 0.5*h*K1);
    K3 = F(t(i) + 0.5*h, y(i) + 0.5*h*K2);
    K4 = F(t(i) + h, y(i) + h*K3);
    T4 = (K1 + 2*K2 + 2*K3 + K4)/6;

    %Saving the next values.
    y(i+1) = y(i) + h*T4;
    t(i+1) = 0 + i*h;
 end

%Plotting th graph
plot(t, y,'Linewidth',2,'DisplayName','Runge Kutta Method')
hold off
legend('location','northwest')
grid on

%}






