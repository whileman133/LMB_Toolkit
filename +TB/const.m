classdef const
    %CONST Constants for the LMB toolbox.

    properties(Constant)
        OCPROOT = TB.getAbsPath('OCP');
        NLEISROOT = TB.getAbsPath('NLEIS');
        MDPROOT = TB.getAbsPath('MDP');
        RPTROOT = TB.getAbsPath('RPT');
        VERSION = '1.0.0';
        AUTHOR = 'Wesley Hileman <whileman@uccs.edu>';

        % Physical constants.
        R = 8.3144598;      % Molar gas constant [J/mol K]
        F = 96485.3329;     % Faraday constant [C/mol]
    end

    methods(Static)
        function f = f(TdegC)
            %F Compute f=F/(R*T) for the supplied T in [C].
            if ~exist('T','var')
                TdegC = 25;
            end
            f = com.const.F/com.const.R/(TdegC+273.15);
        end
    end

end