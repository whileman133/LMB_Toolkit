% function [ROM,tfData] = xraHRA(cellData,xraData,SOC,T)
% 
% Inputs:
%   cellData = cell-parameter data structure loaded with "loadCellParams"
%   xraData  = xra-control parameters loaded with "loadXRA"
%   SOC      = cell state of charge setpoint (between 0 and 1) for the
%              model to be created 
%   T        = cell temperature (in K) for the model to be created
%
% Outputs:
%   ROM      = a structure containing the ROM data
%
% This function creates a reduced-order model for the cell defined by
% "cellData" using the XRA tuning parameters in "xraData" for setpoint
% defined by (SOC,T)... using the hybrid realization algorithm (HRA).
% 
% Copyright (c) 2021 by Gregory L. Plett and M. Scott Trimboli of the
% University of Colorado Colorado Springs (UCCS). This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Intl.
% License, v. 1.0. It is provided "as is", without express or implied
% warranty, for educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L. and Trimboli,
% M. Scott, "Battery Management Systems, Volume III, Physics-Based
% Methods," Artech House, 2021. 

function [ROM,tfData] = xraHRA(cellData,xraData,SOC,T)
  % ----------------- Load HRA tuning parameters -----------------
  Ts   = xraData.Tsamp;    % sampling period [s]
  n1   = xraData.n;       % initial # non-integrator poles
  n    = xraData.nfinal;        % final # non-integrator poles
  kmax = xraData.HRA_Kmax; % number of aliases to consider
  ii   = xraData.HRA_ii;   % self-adjusts...require n <= (ii+1)*numOut

  wnyquist = pi/Ts;          % the Nyquist smapling frequency [rad/s]
  W = unique(xraData.HRA_W)*wnyquist; % Nqyuist rate
  npts = length(W);          % Need to be >> n

  % ----------------- Find discrete-time freq. resp. -----------------
  % Compute continuous-time frequencies needed for general method
  M = (ones(2*kmax+1,1)*W) + 2*pi/Ts*((-kmax:kmax)'*ones(1,length(W)));
  s = 1j*M(:).';

  % Evaluate model parameters at the setpoint
  cellData = evalSetpoint(cellData,s,SOC,T);

  % Evaluate frequency response of all desired TFs
  tflist = xraData.tf;
  [freqResp,hfGain,res0,tfData] = evalTF(tflist,s,cellData);

  % Check if ii is large enough
  numOut = size(freqResp,1);
  if (ii+1)*numOut < n1
    warning('ii not large enough... adjusting');
    ii = ceil(n1/numOut - 1);
  end

  % Remove D term for now; add back into final model
  freqResp = (freqResp - diag(hfGain)*ones(numOut,(2*kmax+1)*npts)).';
  
  % find TFs that are identically zero; delete them for now
  zeroTest = sum(abs(freqResp));
  numTFs = numOut; % number of TFs before deleting any
  indNonzero = find(zeroTest ~= 0);
  indZero = find(zeroTest == 0);
  if ~isempty(indZero)
    freqResp(:,indZero) = [];
    hfGain(indZero) = [];
    res0(indZero) = [];
    numOut = length(hfGain);
  end

  % Set dc gain to 1 (add back into final model)
  dcGain = real(freqResp(M(:)==0,:)); % This is dcgain without D term
  dcGain(dcGain==0) = 1; 
  dcGain(isinf(dcGain)) = 1; 
  dcGain(isnan(dcGain)) = 1;
  freqResp(M(:)==0,:) = dcGain;

  K = diag(dcGain);
  freqResp = freqResp/K;

  % Add together frequency alises to create discrete-time freq resp
  Geval = [];
  for k = 1:numOut
    newData = freqResp(:,k).*(1-exp(-1j*M(:)*Ts))./(1j*M(:));
    newData(M(:)==0) = freqResp(M(:)==0,k)*Ts;
    Geval = [Geval, reshape(newData.',size(M))]; %#ok<AGROW>
  end
  Gdest = (1/Ts)*sum(Geval);
  Gdest = reshape(Gdest.',length(W),numOut).';

  % ------------------------ Execute HRA ------------------------
  % At this point, W is the frequency vector; Gdest is Gd(exp(j*w*Ts)).
  % Need to assemble Y = [Gd; exp(j*w*Ts)*Gd; exp(2j*w*Ts)*Gd; ...]
  % and U = [exp(0j*w*Ts); exp(1j*w*Ts); exp(2j*w*Ts); ...]
  % Then, find [U;G] = [L11 0; L21 L22] * [Q1; Q2]
  % Then, [U1 U2] * [S11 0; 0 0] * [V1'; V2'] = L22
  % Then, Ok = U1*sqrt(S11)
  % Then, C = Ok(1:p,1:n)
  % Finally, [Ok(1:p*(k-1),1:n)] * A = [Ok(p+1:k*p,1:n)].

  % Make data matrices
  ind = 0:ii; 
  ind = repmat(ind,numOut,1);
  U = exp(diag(1j*Ts*(ind(:))')*repmat(W,numOut*(ii+1),1));
  Y = repmat(Gdest,ii+1,1).*U;

  % Perform LQ decomposition of [U;Y]
  LQmat = [real(U) imag(U); real(Y) imag(Y)];
  [~,R2] = qr(LQmat',0); L2 = R2'; w = size(L2);
  L22 = L2(w(1)/2+1:end,w(2)/2+1:end);

  % Find L22 via the LQ decomposition 
  [U,S,~] = svd(L22);
  U1 = U(1:end,1:n1); S11 = S(1:n1,1:n1);
  Ok = real(U1*sqrt(S11));

  % Form A (order n),and D matrices
  A = Ok(1:numOut*ii,1:n1)\Ok(numOut+1:(ii+1)*numOut,1:n1);
  D = hfGain;

  % Diagonalize A and replace unstable poles by their reciprocals
  % and complex poles by their magnitude
  eigA = eig(A);
  eigA(abs(eigA)>1) = 1./eigA(abs(eigA)>1); % flip unstable
  eigA = unique(sort(eigA,'descend'));
  A = diag(abs(eigA)); % make sure real and positive
  diagA = diag(A);

  % ----------------- Compute C matrix (order n) -----------------
  % Use Lagrange multiplier method to compute C to enforce dc gain
  C = zeros(numOut,length(diagA));
  for kk = 1:numOut
    G0 = real(Gdest(kk,find(W==0,1,'first'))); % should be real...
    M1 = zeros(1,length(diagA)); 
    M2 = 0*A;
    for k = 1:length(W) % W,Gdest
      Mx = 1./(exp(1j*W(k)*Ts)-diagA);
      M1 = M1 + real(Gdest(kk,k)*Mx');
      M2 = M2 + real(Mx*Mx');
    end
    M0 = 1./(1-diagA);
    M  = [2*M2' -M0; M0' 0];
    N  = [2*M1'; G0'];
    CL = pinv(M)*N;
    C(kk,:) = real(CL(1:end-1,:))'; % should be real, but sometimes isn't\
  end
  B = ones(length(diagA),1); % order n

  % ----------- Balanced reduction, eliminate complex poles; rescale C -------
  sysHRA = ss(A,B,K*C,D,Ts);
  opt = balredOptions('StateElimMethod','Truncate');
  sysHRA = balred(sysHRA,n,opt);
  sysHRA = canon(sysHRA,'modal',Inf);
  [A,~,~,~] = ssdata(sysHRA);

  % Again, diagonalize A and replace unstable poles by their reciprocals
  % and complex poles by their magnitude
  eigA = eig(A);
  if isfield(xraData,'debug')
    if xraData.debug
      if any(abs(eigA) >= 1), cprintf([1,1/2,0],' - Unstable A... corrected\n'); end
      if any(real(eigA) < 0), cprintf([1,1/2,0],' - Oscillatory A... corrected\n'); end
      if any(eigA ~= conj(eigA)), cprintf([1,1/2,0],' - Complex A... corrected\n'); end
    end
  end
  
  eigA(abs(eigA)>1) = 1./eigA(abs(eigA)>1); % flip unstable
  eigA = abs(sort(eigA,'descend'));
  A = diag(eigA); % make sure real and positive
  diagA = diag(A);

  % ----------------- Compute C matrix again (order nfinal) ---------------
  % Use Lagrange multiplier method to compute C to enforce dc gain
  C = zeros(numOut,n);
  for kk = 1:numOut
    G0 = real(Gdest(kk,find(W==0,1,'first'))); % should be real...
    M1 = zeros(1,n); 
    M2 = 0*A;
    for k = 1:length(W) % W,Gdest
      Mx = 1./(exp(1j*W(k)*Ts)-diagA);
      M1 = M1 + real(Gdest(kk,k)*Mx');
      M2 = M2 + real(Mx*Mx');
    end
    M0 = 1./(1-diagA);
    M = [2*M2' -M0; M0' 0];
    N = [2*M1'; G0'];
    CL = pinv(M)*N;
    C(kk,:) = real(CL(1:end-1,:))'; % should be real, but sometimes isn't
  end
  B = ones(n,1);
  C = K*C;

  % Add back in integrator ...
  if ~prod(res0 == 0)       % Note: res0 is cts-time residue, not disc-time
    A(n+1,n+1) = 1; B = [B;1]; C = [C res0*Ts]; %20211102: Add factor of Ts
  end
  
  % Add back in deleted TFs...
  if ~isempty(indZero)
    Csmall = C; Dsmall = D;
    C = zeros(numTFs,size(Csmall,2));
    D = zeros(numTFs,1);
    C(indNonzero,:) = Csmall;
    D(indNonzero,:) = Dsmall;
  end

  % Save to structure
  ROM.xRA = 'HRA';
  ROM.T   = T;
  ROM.SOC = SOC;
  ROM.A   = A;
  ROM.B   = B;
  ROM.C   = C;
  ROM.D   = D;
  
  if isfield(xraData,'debug')
    if xraData.debug
      % plotting code... make system w/o integrator
      sysHRA = ss(A(1:end-1,1:end-1),B(1:end-1),C(:,1:end-1),D,Ts);

      % Get frequency response of xra-realized system
      [magEst,phaseEst] = bode(sysHRA,W);
      magEst = squeeze(magEst); phaseEst = squeeze(phaseEst);

      % Add dcgain and D term back into true freq response for plotting 
      Gdest = K*Gdest + diag(D)*ones(numOut,npts);
      magTrue = abs(Gdest); phaseTrue = angle(Gdest);

      % Condition phaseEst, phaseTrue so they start in same multiple of 2*pi
      % Note that "unwrap" works on radian angles, hence the conversions...
      phaseTrue = unwrap(phaseTrue,[],2)*180/pi;
      phaseEst = unwrap(phaseEst*pi/180,[],2)*180/pi;
      k = round((phaseEst(:,2)-phaseTrue(:,2))/360);
      phaseEst = phaseEst - k(:,ones(1,size(phaseEst,2)))*360;

      figure; clf;
      subplot(1,2,1);
      semilogx(W,20*log10(magTrue),'linewidth',1.5); hold on
%       semilogx(W,20*log10(magTrue(4,:)),'linewidth',1.5); hold on
      set(gca,'colororderindex',1);       
      semilogx(W,20*log10(abs(magEst)),'x'); hold on;
%       semilogx(W,20*log10(abs(magEst(4,:))),'x'); hold on;
      semilogx([pi/Ts pi/Ts],ylim,'k:');
      xlim([min(W) 1.1*max(W)]); grid on
      xlabel('Analog frequency (rad/s)'); ylabel('Magnitude (dB)');
      title(sprintf('Magnitude plots: Line = truth; x = est.; SOC = %d%%',100*SOC));

      subplot(1,2,2);
      semilogx(W,phaseTrue,'linewidth',1.5); hold on
%       semilogx(W,phaseTrue(4,:),'linewidth',1.5); hold on
      set(gca,'colororderindex',1);       
      semilogx(W,phaseEst,'x'); hold on;
%        semilogx(W,phaseEst(4,:),'x'); hold on;
      semilogx([pi/Ts pi/Ts],ylim,'k--');
      xlim([min(W) 1.1*max(W)]); grid on
      xlabel('Analog frequency (rad/s)'); ylabel('Phase (deg)');
      title(sprintf('Phase plots: Line = truth; x = est.; SOC = %d%%',100*SOC));
      drawnow

      save = sprintf('figures/debugging/debugging_psiT_%d',100*SOC);
      saveas(gcf,save,'png')

    end
  end
end

