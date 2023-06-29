classdef DiagonalEstimate < handle
    %DIAGONALOCPESTIMATE Implements the diagonal interpolation method for
    %  estimating OCP.
    %
    % To convert to a regular OcpEstimate object, call one of the
    % following:
    %   useAvg() to use the OCP estimate obtained by averaging shifted dis/charge curves.
    %   useDis() to use the OCP estimate obtained by shifting discharge curve.
    %   useChg() to use the OCP estimate obtained by shifting charge curve.

    properties(SetAccess=protected)
        avgV   % OCP estimate obtained by averaging shifted dis/charge curves.
        disV   % OCP estimate obtained by shifting discharge curve.
        chgV   % OCP estimate obtained by shifting charge curve.
        comZ   % Common stiochiometry vector for all estimates.

        disdzdvUnifV   % Differential capacity, discharge curve, over uniform voltage vector.
        chgdzdvUnifV   % Differential capacity, charge curve, over uniform voltage vector.
        disdzdvUnifZ   % Differential capacity, discharge curve, over uniform stiochiometry vector.
        chgdzdvUnifZ   % Differential capacity, charge curve, over uniform stiochiometry vector.
        unifV          % Uniform potential vector for differential capacity curves.
        unifZ          % Uniform stiochiometry vector for differential capacity curves.

        lagZ      % Amount by which dis/charge curves shifted over stiochiometry [unitless].
        lagV      % Amount by which dis/charge curves shifted over potential [V].
 
        ocptest   % OCP test object from which estimate is derived.
        dvbin     % Size of OCP bins used in differential capacity calculation [V].
    end

    methods
        function obj = DiagonalEstimate(ocptest, dvbin, varargin)
            %DIAGONALOCPESTIMATE Estimate OCP from the two dis/charge
            %  curves using the diagonal interpolation method.
            %
            % ocpTest = OcpTest object containing the dis/charge data.
            % dvbin = size of OCP bins to use in differential capacity
            %   computation (histogram method.)

            p = inputParser;
            p.addOptional('minV',[]);
            p.addOptional('maxV',[]);
            p.addOptional('maxLagV',0.1);
            p.addOptional('minZ',[],@(x)0<=x&&x<=1);
            p.addOptional('maxZ',[],@(x)0<=x&&x<=1);
            p.addOptional('maxLagZ',0.1,@(x)0<=x&&x<=1);
            p.parse(varargin{:});
            minV = p.Results.minV;
            maxV = p.Results.maxV;
            maxLagV = p.Results.maxLagV;
            minZ = p.Results.minZ;
            maxZ = p.Results.maxZ;
            maxLagZ = p.Results.maxLagZ;
            
            nsamples = ceil((length(ocptest.disZ)+length(ocptest.chgZ))/2);
            comZ = linspace(0, 1, nsamples)';

            % Interpolate dis/charge curves over common relative 
            % stiochiometry vector, linearly extrapulate to fill entire 
            % 0 to 1 range. 
            disV = interp1(ocptest.disZ, ocptest.disV, comZ, 'linear', 'extrap');
            chgV = interp1(ocptest.chgZ, ocptest.chgV, comZ, 'linear', 'extrap');

            % Compute differential capacity.
            [disdzdv_V, disVout, disdzdv_Z, disZout] = ocp.smoothdiff(comZ, disV, dvbin);
            [chgdzdv_V, chgVout, chgdzdv_Z, chgZout] = ocp.smoothdiff(comZ, chgV, dvbin);

            % Interpolate differential capacity over uniform stiochiometry
            % vector (required to compute cross-correlation w/r/t stiochiometry below.)
            if isempty(minZ), minZ = max(min(disZout), min(chgZout)); end
            if isempty(maxZ), maxZ = min(max(disZout), max(chgZout)); end
            nsamples = ceil((length(disZout)+length(chgZout))/2);
            unifZ = linspace(minZ, maxZ, nsamples)';
            unifdz = mean(diff(unifZ));
            disdzdvUnifZ = interp1(disZout, disdzdv_Z, unifZ, 'linear', 'extrap');
            chgdzdvUnifZ = interp1(chgZout, chgdzdv_Z, unifZ, 'linear', 'extrap');

            % Interpolate differential capacity over uniform potential
            % vector (required to compute cross-correlation w/r/t potential below.)
            if isempty(minV), minV = max(min(disVout), min(chgVout)); end
            if isempty(maxV), maxV = min(max(disVout), max(chgVout)); end
            nsamples = ceil((length(disVout)+length(chgVout))/2);
            unifV = linspace(minV, maxV, nsamples)';
            unifdv = mean(diff(unifV));
            disdzdvUnifV = interp1(disVout, disdzdv_V, unifV, 'linear', 'extrap');
            chgdzdvUnifV = interp1(chgVout, chgdzdv_V, unifV, 'linear', 'extrap');

            % Compute cross-correlation and lag w/r/t stiochiometry.
            maxLag = ceil(maxLagZ/unifdz);
            [correlationZ, lagZ] = xcorr(disdzdvUnifZ, chgdzdvUnifZ, maxLag);
            % ! Normalize to account for finite length of the vectors.
            correlationZ = correlationZ./(length(disdzdvUnifZ)-abs(lagZ))';
            [~, idxMax] = max(correlationZ);
            lagZ = abs(lagZ(idxMax)*unifdz);

            % Compute cross-correlation and lag w/r/t potential.
            maxLag = ceil(maxLagV/unifdv);
            [correlationV, lagV] = xcorr(disdzdvUnifV, chgdzdvUnifV, maxLag);
            % ! Normalize to account for finite length of the vectors.
            correlationV = correlationV./(length(disdzdvUnifV)-abs(lagV))';
            [~, idxMax] = max(correlationV);
            lagV = abs(lagV(idxMax)*unifdv);

            % Perform the diagonal interpolation by shifting dis/charge
            % curves according to the observed lags.
            disShiftZ = ocptest.disZ + lagZ/2;
            disShiftV = ocptest.disV + lagV/2;
            chgShiftZ = ocptest.chgZ - lagZ/2;
            chgShiftV = ocptest.chgV - lagV/2;

            % Interpolate over common stiochiometry vector.
            disInterpV = interp1(disShiftZ, disShiftV, comZ, 'linear', 'extrap');
            chgInterpV = interp1(chgShiftZ, chgShiftV, comZ, 'linear', 'extrap');

            % Average shifted dis/charge curves.
            interpV = (disInterpV + chgInterpV)/2;

            % Store results into object.
            obj.ocptest = ocptest;
            obj.avgV = interpV;
            obj.disV = disInterpV;
            obj.chgV = chgInterpV;
            obj.comZ = comZ;
            obj.disdzdvUnifV = disdzdvUnifV;
            obj.chgdzdvUnifV = chgdzdvUnifV;
            obj.disdzdvUnifZ = disdzdvUnifZ;
            obj.chgdzdvUnifZ = chgdzdvUnifZ;
            obj.unifV = unifV;
            obj.unifZ = unifZ;
            obj.lagV = lagV;
            obj.lagZ = lagZ;
            obj.dvbin = dvbin;
        end

        function avgEst = useAvg(obj)
            %USEAVG Build an OcpEstimate using the OCP estimate obtained by
            %  averaging shifted dis/charge curves.
            avgEst = ocp.Estimate(obj.ocptest,obj.dvbin,obj.comZ,obj.avgV,'Diagonal, Averaged');
        end

        function disEst = useDis(obj)
            %USEDIS Build an OcpEstimate using the OCP estimate obtained by
            %  shifting the discharge curve.
            disEst = ocp.Estimate(obj.ocptest,obj.dvbin,obj.comZ,obj.disV,'Diagonal, Shifted Discharge');
        end

        function chgEst = useChg(obj)
            %USECHG Build an OcpEstimate using the OCP estimate obtained by
            %  shifting the charge curve.
            chgEst = ocp.Estimate(obj.ocptest,obj.dvbin,obj.comZ,obj.chgV,'Diagonal, Shifted Charge');
        end

        function figures = plot(obj)
            %PLOT Plot estimate results.

            f1 = figure;
            plot(obj.comZ, obj.avgV); hold on;
            plot(obj.ocptest.disZ, obj.ocptest.disV, 'k');
            plot(obj.ocptest.chgZ, obj.ocptest.chgV, 'k');
            ylim([obj.ocptest.vmin obj.ocptest.vmax]);
            title(sprintf('%s OCP Estimate, %.2f\\circC', obj.ocptest.name, obj.ocptest.temp));
            xlabel('Relative Lithiation');
            ylabel('Potential versus Li/Li^+ [V]');
            legend('Diagonal Estimate', 'Dis/charge Data');
            thesisFormat([0.1 0.05 0.05 0]);

            f2 = figure;
            plot(obj.comZ, obj.disV); hold on;
            plot(obj.comZ, obj.chgV);
            plot(obj.ocptest.disZ, obj.ocptest.disV, 'k');
            plot(obj.ocptest.chgZ, obj.ocptest.chgV, 'k');
            ylim([obj.ocptest.vmin obj.ocptest.vmax]);
            title(sprintf('%s Dis/charge OCP Estimates, %.2f\\circC', obj.ocptest.name, obj.ocptest.temp));
            xlabel('Relative Lithiation');
            ylabel('Potential versus Li/Li^+ [V]');
            legend('Discharge', 'Charge');
            thesisFormat([0.1 0.05 0.05 0]);

            f3 = figure;
            plot(obj.unifV, obj.disdzdvUnifV); hold on;
            plot(obj.unifV, obj.chgdzdvUnifV);
            xlim([obj.ocptest.vmin obj.ocptest.vmax]);
            title(sprintf('%s Diff. Capacity vs OCP, %.2f\\circC', obj.ocptest.name, obj.ocptest.temp));
            xlabel('Potential versus Li/Li^+ [V]');
            ylabel('Differential Capacity [V^{-1}]')
            legend('Discharge', 'Charge');
            thesisFormat([0.1 0.05 0.05 0]);

            f4 = figure;
            plot(obj.unifZ, obj.disdzdvUnifZ); hold on;
            plot(obj.unifZ, obj.chgdzdvUnifZ);
            title(sprintf('%s Diff. Capacity vs Lithiation, %.2f\\circC', obj.ocptest.name, obj.ocptest.temp));
            xlabel('Relative Lithiation');
            ylabel('Differential Capacity [V^{-1}]')
            legend('Discharge', 'Charge');
            thesisFormat([0.1 0.05 0.05 0]);

            figures = [f1 f2 f3 f4];
        end
    end
end