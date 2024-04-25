% compareNMC811Ds.m
clear; close all; clc;
addpath('..');
TB.addpaths;

RuessDs = readtable('NMC811-Ds-RuessEtAl.xlsx');
UocpREF = RuessDs.Uocp;
DsREF_cm2secN1 = RuessDs.Ds;

figure;
semilogy(UocpREF,DsREF_cm2secN1,'o-');
thesisFormat;