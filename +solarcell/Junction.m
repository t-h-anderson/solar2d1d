classdef Junction
% solarcell.Junction
%
% Value class wrapping the configuration produced by DesignJunction.m.
% This is the first step of the OOP migration; it does not change any
% downstream behaviour. DesignJunction.m stays in place as the
% regression oracle, and the snapshot test in tests/JunctionSnapshotTest
% asserts that solarcell.Junction().toStruct() is isequaln-equal to it.
%
% Once the snapshot test is locked in we can move callers off
% DesignJunction one at a time and start adding real behaviour to this
% class (parameter sweeps, validation, named fixtures, ...).
%
% Layout of the output struct (matches DesignJunction):
%   .input             nmLx, VVext, icm3isG
%   .material          stack description (nsec sections)
%   .optical.periodic  grating description
%   .phys              physical constants + derived scaling factors

    properties
        % --- Sim-time inputs (sim.input) ---
        nmLx        (1,1) double = 500          % Grating period (nm)
        VVext       (1,1) double = 0.0          % External voltage (V)
        icm3isG     (1,1) double = -inf         % Generation rate; -inf = take from RCWA

        % --- Material stack (sim.material) ---
        nsec        (1,1) double = 9
        elecsec     (1,:) double = [0,0,0,1,1,1,0,0,0]
        periodicsec (1,:) double = [0,0,0,0,0,0,0,1,0]
        nmSec       (1,:) double = [100,110,100,10,25,20,20,1,100]
        icm3ND      (1,:) double = [0,0,0,1e17,5e17,-1e16,0,0,0]
        material    char         = 'Air & MgF2 & AZO &ZnO& CdS & CIGS & AZO &AZO& Silver'
        VEg0        (1,:) double = [0,0,0,3.3,2.4,1.13,0,0,0]
        epsdc       (1,1) double = 11.68
        icm3isalpha (1,1) double = 1e-10
        istausrhn   (1,1) double = 0e-12
        istausrhp   (1,1) double = 0e-9
        Cn          (1,1) double = 0e-20
        Cp          (1,1) double = 0e-20
        ET          (1,1) double = 0.5
        eps         (1,:) double = [1, 1, 1, 1, 1, 15+2i, 1, 1, 1]

        % --- Optical grating (sim.optical.periodic) ---
        zr_flag     (1,:) double = 1
        zeta        (1,1) double = 0.5
        relief      char         = '@(zeta) cos(2*pi*zeta)'
        material2   char         = 'Silver'
        eps2        (1,1) double = -10.4356+0.8294i

        % --- Tunable scaling input ---
        % DesignJunction hard-codes this at 1.5; exposing it as a property
        % so future callers can override without editing the class.
        cm2iVismus  (1,1) double = 1.5
    end

    methods
        function sim = toStruct(obj)
            % Materialise the struct expected by every legacy consumer
            % (DesignSim, BuildOptInput, RunOptical, RunElectrical, ...).
            % Field layout matches DesignJunction.m exactly.

            obj.validate();

            sim.input    = obj.inputStruct();
            sim.material = obj.materialStruct();
            sim.optical.periodic = obj.gratingStruct();
            sim.phys     = obj.physStruct();
        end
    end

    methods (Access = private)
        function validate(obj)
            if length(obj.nmSec) ~= obj.nsec
                error('solarcell:Junction:lengthMismatch', ...
                    'Lengths not defined for each section of junction');
            end
        end

        function s = inputStruct(obj)
            s = struct( ...
                'nmLx',    obj.nmLx, ...
                'VVext',   obj.VVext, ...
                'icm3isG', obj.icm3isG);
        end

        function s = materialStruct(obj)
            s = struct( ...
                'nsec',        obj.nsec, ...
                'elecsec',     obj.elecsec, ...
                'periodicsec', obj.periodicsec, ...
                'nmSec',       obj.nmSec, ...
                'icm3ND',      obj.icm3ND, ...
                'material',    obj.material, ...
                'VEg0',        obj.VEg0, ...
                'epsdc',       obj.epsdc, ...
                'icm3isalpha', obj.icm3isalpha, ...
                'istausrhn',   obj.istausrhn, ...
                'istausrhp',   obj.istausrhp, ...
                'Cn',          obj.Cn, ...
                'Cp',          obj.Cp, ...
                'ET',          obj.ET, ...
                'eps',         obj.eps);
        end

        function s = gratingStruct(obj)
            s = struct( ...
                'zr_flag',   obj.zr_flag, ...
                'zeta',      obj.zeta, ...
                'relief',    obj.relief, ...
                'material2', obj.material2, ...
                'eps2',      obj.eps2);
        end

        function phys = physStruct(obj)
            % Mirrors the order of assignments in DesignJunction.m so the
            % field-name set is identical. (isequaln on structs ignores
            % field order, but matching it keeps diffs minimal during the
            % migration.)
            phys.Cq           = 1.6021766e-19;
            phys.Jshbar       = 1.0545718e-34;
            phys.VVth         = 300 * 8.61733034e-5;
            phys.CiVimeps0    = 8.85418782e-12;
            phys.CiVicmeps0   = phys.CiVimeps0 / 100;
            phys.mkgis2iA2mu0 = 4*pi*10^(-7);
            phys.Oeta0        = 376.73031346177;
            phys.misc         = 299792458;
            phys.m2kgish      = 6.6260696e-34;
            phys.Enormconst   = sqrt(1 / (phys.CiVimeps0 * phys.misc));

            phys.nmLs         = obj.nmSec * obj.elecsec.';
            phys.cmLs         = cm_from_nm(phys.nmLs);
            phys.icm3Ns       = max([max(abs(obj.icm3ND)), 1e14]);
            phys.VVs          = phys.VVth;
            phys.cm2iVismus   = obj.cm2iVismus;
            phys.icm3isGs     = phys.cm2iVismus * phys.VVs * phys.icm3Ns / phys.cmLs^2;

            phys.cm3isalphas  = (phys.cm2iVismus * phys.VVs) / (phys.cmLs^2 * phys.icm3Ns);
            phys.isTaus       = phys.cmLs^2 / (phys.cm2iVismus * phys.VVs);
            phys.cm6isCs      = (phys.cm2iVismus * phys.VVs) / (phys.cmLs^2 * phys.icm3Ns^2);

            % DesignJunction reassigns cmLs here. Reproduce verbatim so a
            % future caller comparing field counts gets the same answer.
            phys.cmLs         = cm_from_nm(phys.nmLs);
            phys.cm2isDs      = phys.cm2iVismus * phys.VVs;
            phys.mAicm2Js     = 1000 * phys.Cq * phys.VVs * phys.cm2iVismus * phys.icm3Ns / phys.cmLs;
            phys.icm2isGRs    = phys.VVs * phys.cm2iVismus * phys.icm3Ns / phys.cmLs^2;
            phys.lambda2s     = phys.CiVicmeps0 * phys.VVth / (phys.cmLs^2 * phys.Cq * phys.icm3Ns);
            phys.Qscale       = phys.Oeta0 * phys.CiVimeps0 / (phys.Jshbar * phys.Enormconst^2);
        end
    end
end
