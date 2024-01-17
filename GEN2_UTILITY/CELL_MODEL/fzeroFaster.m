% -------------------------------------------------------------------------
% Replace MATLAB's built-in FZERO command with this stripped-down version
% that removes a lot of error checking and tracing capabilities in
% preference for speed. Runs about twice as fast. See help on FZERO
% -------------------------------------------------------------------------
function b = fzeroFaster(FunFcn,x,options,varargin) %#ok<INUSL>
  % I always call with four inputs: anon function, initial guess, empty
  % options [], and a second argument to the anonymous function.
  % I always call with one output, so other outputs from FZERO removed.

  % Initialization
  fcount = 0;
  iter = 0;
  intervaliter = 0;
  tol = eps;

  % Put first feval in try catch
  try
    fx = FunFcn(x,varargin{:});
  catch ME
    if ~isempty(Ffcnstr)
      error('MATLAB:fzero:InvalidFunctionSupplied',...
          getString(message('MATLAB:optimfun:fzero:InvalidFunctionSupplied',...
          sprintf('%s ==> %s','function_handle',Ffcnstr),ME.message)));
    else
      error('MATLAB:fzero:InvalidFunctionSupplied',...
          getString(message('MATLAB:optimfun:fzero:InvalidFunctionSupplied',...
          'function_handle',ME.message)));
    end
  end
  fcount = fcount + 1;
  if fx == 0
    b = x;
    return
  elseif ~isfinite(fx) || ~isreal(fx)
    error('MATLAB:fzero:ValueAtInitGuessComplexOrNotFinite',...
        getString(message('MATLAB:optimfun:fzero:ValueAtInitGuessComplexOrNotFinite')));
  end

  if x ~= 0
    dx = x/50;
  else
    dx = 1/50;
  end

  % Find changes of sign.
  twosqrt = sqrt(2);
  a = x; fa = fx; b = x; fb = fx;

  while (fa > 0) == (fb > 0)
    intervaliter = intervaliter + 1;
    dx = twosqrt*dx;
    a = x - dx;  fa = FunFcn(a,varargin{:});
    fcount = fcount + 1;
    if ~isfinite(fa) || ~isreal(fa) || ~isfinite(a)
        b = NaN; return
    end

    if (fa > 0) ~= (fb > 0) % check for different sign
      break
    end

    b = x + dx;  fb = FunFcn(b,varargin{:});
    if ~isfinite(fb) || ~isreal(fb) || ~isfinite(b)
      b = NaN; return
    end
    fcount = fcount + 1;
  end % while

  fc = fb;
  while fb ~= 0 && a ~= b
    % Ensure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
      c = a;  fc = fa;
      d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
      a = b;    b = c;    c = a;
      fa = fb;  fb = fc;  fc = fa;
    end

    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(m) <= toler) || (fb == 0.0)
      break
    end

    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
      % Bisection
      d = m;  e = m;
    else
      % Interpolation
      s = fb/fa;
      if (a == c)
        % Linear interpolation
        p = 2.0*m*s;
        q = 1.0 - s;
      else
        % Inverse quadratic interpolation
        q = fa/fc;
        r = fb/fc;
        p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
        q = (q - 1.0)*(r - 1.0)*(s - 1.0);
      end
      if p > 0
        q = -q;
      else
        p = -p;
      end
      % Is interpolated point acceptable
      if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
        e = d;  d = p/q;
      else
        d = m;  e = m;
      end
    end % Interpolation

    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler
      b = b + d;
    elseif b > c
      b = b - toler;
    else
      b = b + toler;
    end
    fb = FunFcn(b,varargin{:});
    fcount = fcount + 1;
    iter = iter + 1;
  end % Main loop
end