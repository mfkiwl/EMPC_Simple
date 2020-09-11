function uuout = MPCController(currentx,currentV,currentPg,currentPw,t,currentA,currentB)
global Parameters
persistent Controller
Ad=1;
Bd=Parameters.Ts*[1 -1]*Parameters.ScaleP/Parameters.ScaleK*Parameters.G^2; %[Pw; Pg]
Cd=[1];
Dd=[0];
Np  = 90;
if t==0
    currentx=Parameters.J/2*Parameters.Omegagmin^2*1.1;
    currentV=5.3;
    currentPgg=0.5e6/Parameters.ScaleP;
    currentPww=2.5e6/Parameters.ScaleP;
currentA=zeros(Parameters.Nconst,1);
currentB=ones(Parameters.Nconst,1)*10;
% Define data for MPC controller

% Avoid explosion of internally defined variables in YALMIP
yalmip('clear')

% Setup the optimization problem
u = sdpvar(repmat(2,1,Np),repmat(1,1,Np));
x = sdpvar(repmat(1,1,Np+1),repmat(1,1,Np+1));
sdpvar v
AA = sdpvar(Parameters.Nconst,1);
BB = sdpvar(Parameters.Nconst,1);
sdpvar PG
sdpvar PW
%% find windspeed bin and fraction
% qq=find(sign(diff(sign(Parameters.BOUNDPW.V-currentV))),1);
% fraction=1-(currentV-Parameters.BOUNDPW.V(qq))/(Parameters.BOUNDPW.V(qq+1)-Parameters.BOUNDPW.V(qq));

%% Setup Lp+QP
constraints = [];
objective = 0;
tic
currentxx=currentx/Parameters.ScaleK;
currentPgg=currentPg/Parameters.ScaleP;
currentPww=currentPw/Parameters.ScaleP;
for k = 1:Np
    objective = objective - u{k}(2); % optimize the available power
    
    constraints = [constraints, x{k+1}== Ad*x{k}+Bd*u{k}];
    constraints = [constraints, 0 <= u{k}(2)<= sqrt(2*x{k+1}*Parameters.ScaleK/Parameters.J)*Parameters.Tmax/Parameters.ScaleP];
    constraints = [constraints, 0 <= u{k}(2)<= 5e6/Parameters.ScaleP];
    
    constraints = [constraints,  Parameters.J/2*Parameters.Omegagmin^2/Parameters.ScaleK<= (x{k+1})<= Parameters.J/2*Parameters.Omegagmax^2/Parameters.ScaleK];
    for ijj=1:1:Parameters.Nconst
        constraints = [constraints, 0<= (u{k}(1))<= (AA(ijj)+BB(ijj)*x{k+1})*v^3];
    end
    
    if k==1
        constraints = [constraints, -0.2e6*Parameters.Ts/Parameters.ScaleP <= (u{k}(2)-PG)<= 0.2e6*Parameters.Ts/Parameters.ScaleP];
    else
        constraints = [constraints, -0.2e6*Parameters.Ts/Parameters.ScaleP <= (u{k}(2)-u{k-1}(2))<= 0.2e6*Parameters.Ts/Parameters.ScaleP];
    end
    
        if k==1
        constraints = [constraints, -0.5e6*Parameters.Ts/Parameters.ScaleP <= (u{k}(1)-PW)<= 0.5e6*Parameters.Ts/Parameters.ScaleP];
    else
        constraints = [constraints, -0.5e6*Parameters.Ts/Parameters.ScaleP <= (u{k}(1)-u{k-1}(1))<=0.5e6*Parameters.Ts/Parameters.ScaleP];
    end
end

% Define an optimizer object which solves the problem for a particular
% initial state and reference

%ops = sdpsettings('solver','Mosek','verbose',2);
ops = sdpsettings('solver','Mosek');
Controller = optimizer(constraints,objective,ops,{x{1},v,AA,BB,PG,PW},{[u{:}],[x{:}]});
toc
% And use it here too
[uout_complete] = Controller{{currentxx,currentV,currentA,currentB,currentPgg,currentPww}};
toc

uout=uout_complete{1}(:,1)*Parameters.ScaleP;
solutionP=uout_complete{1}*Parameters.ScaleP; 
solutionX=uout_complete{2}*Parameters.ScaleK; 
%% define outputs
Pw=uout(1);
Pg=uout(2);
omega_g=sqrt(uout_complete{2}(:,1)*Parameters.ScaleK*2/Parameters.J);
Tg=uout(2)/omega_g/Parameters.Effic;

%% compute pitch
pitch=-3:0.01:20;
TSR=omega_g/Parameters.G*Parameters.R/currentV;
Cpp=Pw/(1/2*Parameters.rho*pi*Parameters.R^2*currentV^3)
Cp = max(interp2(Parameters.Rotor_Lamda,Parameters.Rotor_Pitch,Parameters.Rotor_cP,TSR,pitch,'cubic'),0);
[MinObje,IndexPitch]=min((Cp-Cpp).^2)
Pitch=pitch(IndexPitch)
uuout=[Pg; Tg;Pw;Pitch;TSR;Cpp;omega_g];
else
    

currentxx=currentx/Parameters.ScaleK;
currentPgg=currentPg/Parameters.ScaleP;
currentPww=currentPw/Parameters.ScaleP;
tic
% And use it here too
[uout_complete] = Controller{{currentxx,currentV,currentA,currentB,currentPgg,currentPww}};
toc

uout=uout_complete{1}(:,1)*Parameters.ScaleP;
solutionP=uout_complete{1}*Parameters.ScaleP; 
solutionX=uout_complete{2}*Parameters.ScaleK; 
%% define outputs
Pw=uout(1);
Pg=uout(2);
omega_g=sqrt(uout_complete{2}(:,1)*Parameters.ScaleK*2/Parameters.J);
Tg=uout(2)/omega_g/Parameters.Effic;

%% compute pitch
pitch=-3:0.005:20;
TSR=omega_g/Parameters.G*Parameters.R/currentV;
Cpp=Pw/(1/2*Parameters.rho*pi*Parameters.R^2*currentV^3);
Cp = max(interp2(Parameters.Rotor_Lamda,Parameters.Rotor_Pitch,Parameters.Rotor_cP,TSR,pitch,'cubic'),0);
[MinObje,IndexPitch]=min((Cp-Cpp).^2);
Pitch=pitch(IndexPitch);
uuout=[Pg; Tg;Pw;Pitch;TSR;Cpp;omega_g];
end
