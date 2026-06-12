%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%% Test avec le code d'Az*
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
clc
clear
format long


%Parameters

Param.E = 10^4;                       % Young modulus
Param.nu = 0.38;                      % Poisson ratio
Param.L = 0.5;                          % Rod length
Param.r = 0.03;                       % rod radius
Param.A =  pi*Param.r^2;              % Cross section area
Param.V = Param.L*Param.A;            % rod volume
Param.G = Param.E/(2*(1+Param.nu));   % shear modulus
Param.m = 0.05;                        % Beam mass
Param.rho = Param.m/Param.V;          % Mass density
Param.g =[0 ; 0 ; 0 ; 0 ; 0; -9.81]; % gravity
Param.J1 = pi*Param.r^4/2; % Polar inertia moment
Param.J2 = pi*Param.r^4/4; % Inertia moment
Param.J3 = pi*Param.r^4/4; % Inertia moment
Param.H = diag([Param.G*Param.J1 , Param.E*Param.J2 , Param.E*Param.J3 , Param.E*Param.A , Param.G*Param.A , Param.G*Param.A]); %Hooke tensor
Param.M = diag([Param.rho*Param.J1,Param.rho*Param.J2,Param.rho*Param.J3,Param.rho*Param.A,Param.rho*Param.A,Param.rho*Param.A]);%cross-sectional inertia

Param.Forces_Tendons=[0;0;0;0];           % Pas d'actionnement
Param.Ftip= [0;0;0;0;0;0];

Param.dX=Param.L/200;                         % spatial step
Param.n_X=length(Param.dX:Param.dX:Param.L); % Number of section

Param.DeltaX2 = Param.L/200;
Param.n_seg = length(Param.DeltaX2:Param.DeltaX2:Param.L);

Param.Rb = 0.02;     % Distance between a tendon and the backbone
Param.Tendon_coordinate = @(theta, distance) distance*[0; cos(theta); sin(theta)]; % coordinate of the tendon in a cross-section
Param.Tendons_list=zeros(3*(Param.n_seg+1),4);

for i=1:Param.n_seg+1
  Param.Tendons_list(3*(i-1)+1:3*i, :) = [Param.Tendon_coordinate(0, Param.Rb),Param.Tendon_coordinate(pi/2, Param.Rb),Param.Tendon_coordinate(pi,Param.Rb),Param.Tendon_coordinate(3*pi/2, Param.Rb)];

end

Param.n=3;  %strain modes number
Param.na=2; %number of actuated strains

% BC (fixed base)
Param.g0 = eye(4);

% Projection matrices
Param.B     = [0 0;1 0;0 0;0 1;0 0;0 0];
Param.B_bar = [1 0 0 0;0 0 0 0;0 1 0 0;0 0 0 0;0 0 1 0;0 0 0 1];
Param.xi_0=[0;0;0;1;0;0];
Param.xi_a0=Param.B'*Param.xi_0;
Param.xi_c=[0;0;0;0];

Param.Ha=Param.B'*Param.H*Param.B;%the matrix of the reduced Hooke coefficients

%%%% ?????
Y = 0:Param.DeltaX2:Param.L;
SK = zeros(Param.n*Param.na, Param.n*Param.na);
for k = 1:Param.n_seg+1
    phival = Phi(Param.na, Param.n, Y(k), Param.L);
    fK = phival'*Param.Ha*phival;
    if k == 1
        fK0 = fK;
    elseif k == (Param.n_seg+1)
        fKn = fK;
    else
        SK = SK+fK;
    end
end
Param.Keps = Param.DeltaX2.*((fK0+fKn)./2 + SK);

%% Resolution
init = zeros(Param.na*Param.n, 1);

##options = optimset('fsolve');
##options.MaxFunEvals = 1000000000000;
##options.MaxIter = 50000000000;
##options.TolFun=1e-20;
##options.Jacobian='off';
##options.Display = 'iter';

options = optimset('Display', 'iter',
                   'MaxFunEvals', 1e12,
                   'MaxIter', 5e10,
                   'TolFun', 1e-20);

tic
[q1, fval, exitflag] = fsolve(@(var) Static_V2(var, Param), init, options);
toc

final_error = norm(fval)
fprintf('residual norm : %e\n', final_error);


%% Beam reconstruction

r1=zeros(3,Param.n_X+1);
Q1=zeros(4,Param.n_X+1);

Q1(:,1)=[1;0;0;0];

for j=1:Param.n_X
  %   g=Geometric_model((j-1)*Param.dX,q,Param);
  %   r(:,j)=g(1:3,4);
        phi=Phi(Param.na,Param.n,(j-1)*Param.dX,Param.L);%Functions basis values at X
        xia1=Param.xi_a0 + phi*q1;
        xi1=Param.B*xia1+Param.B_bar*Param.xi_c;
        %xi(4:6)=[1;0;0];
        K1=xi1(1:3); %angular strain
        Gamma1=xi1(4:6);%Linear strain
        R1=eye(3) + 2/(Q1(:,j)'*Q1(:,j)) * [-Q1(3,j)^2-Q1(4,j)^2, Q1(2,j)*Q1(3,j)-Q1(4,j)*Q1(1,j),Q1(2,j)*Q1(4,j) + Q1(3,j)*Q1(1,j) ;
                                 Q1(2,j)*Q1(3,j)+Q1(4,j)*Q1(1,j), -Q1(2,j)^2-Q1(4,j)^2,Q1(3,j)*Q1(4,j) - Q1(2,j)*Q1(1,j) ;
                                 Q1(2,j)*Q1(4,j)-Q1(3,j)*Q1(1,j), Q1(3,j)*Q1(4,j) + Q1(2,j)*Q1(1,j), -Q1(2,j)^2-Q1(3,j)^2];
        Q_X1 = [ 0, -K1(1), -K1(2), -K1(3);
                 K1(1), 0, K1(3), -K1(2);
                K1(2), -K1(3), 0, K1(1);
                K1(3), K1(2), -K1(1), 0 ] * Q1(:,j)/2;
        r_X1 = R1*Gamma1;
        Q1(:,j+1)=Q1(:,j) + Q_X1*Param.dX;
        r1(:,j+1)=r1(:,j) + r_X1*Param.dX;
        R1=eye(3) + 2/(Q1(:,j+1)'*Q1(:,j+1)) * [-Q1(3,j+1)^2-Q1(4,j+1)^2, Q1(2,j+1)*Q1(3,j+1)-Q1(4,j+1)*Q1(1,j+1),Q1(2,j+1)*Q1(4,j+1) + Q1(3,j+1)*Q1(1,j+1) ;
                                 Q1(2,j+1)*Q1(3,j+1)+Q1(4,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(4,j+1)^2,Q1(3,j+1)*Q1(4,j+1) - Q1(2,j+1)*Q1(1,j+1) ;
                                 Q1(2,j+1)*Q1(4,j+1)-Q1(3,j+1)*Q1(1,j+1), Q1(3,j+1)*Q1(4,j+1) + Q1(2,j+1)*Q1(1,j+1), -Q1(2,j+1)^2-Q1(3,j+1)^2];
end
    %hold on
    plot(r1(1,:),r1(3,:),'r','LineWidth', 4); title('Cosserat rod');axis([-Param.L 2*Param.L -Param.L Param.L -Param.L Param.L]);grid on; daspect([1 1 1])
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); drawnow
    %hold off
% end

export_pos = r1'
save('~/Beam/beam_pos_matlab.txt', 'export_pos', '-ascii');
disp('End of the computation!! \n');
