%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%%    Description

%     This code provides the fuison results of PANchromatic (PAN) and MultiSpectral (MS) images
%     based on Particle Swarm Optimization (PSO) algorithm proposed in the following reference:

%     [1] A. Azarang and H. Ghassemian, "An adaptive multispectral image fusion
%     using particle swarm optimization," in Proc. Iranian Conf. Elec. Eng.
%     (ICEE), May 2017, pp. 1708-1712.

%     [2] S. Rahmani, M. Strait, D. Merkurjev, M. Moeller, and T. Wittman,
%     "An adaptive IHS Pan-sharpening method," IEEE Geosci. Remote Sens.
%     Lett., vol. 7, no. 4, pp. 746-750, Oct. 2010.

%     [3] Y. Leung, J. Liu, and J. Zhang, "An improved adaptive intensity-huesaturation
%     method for the fusion of remote sensing images," IEEE Geosci.
%     Remote Sens. Lett., vol. 11, no. 5, pp. 985-989, May 2014.
%%    The steps invovled in the proposed method is summerized as:

%     1) Pre-processing the datasets,
%     2) Estimating the primitive detail map,
%     3) Extracting the primitive PAN and MS edge detectors,
%     4) Optimizing the weights of edge detectors using the PSO algorithm
%     and ERGAS loss function to be applied on primitive detail map.

% 全色和多光谱图像融合（Panchromatic and multi-spectral image fusion，Pan sharping）
% 是将源图像的空间和光谱信息融合成一幅融合图像，该融合图像具有比任何源图像更高的空间和光谱分辨率，
% 对下游任务更可靠。它已广泛应用于各种应用的图像解释和预处理。
% 为了通过考虑全色和多光谱图像之间的空间和光谱关系来获得更好的融合结果，
% 已经提出了大量的方法

% 需引用gramm的论文
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading the images and pre-processing steps


addpath Algorithm;
addpath Objective_Evaluation;
addpath Utils;

%% Make the PAN and MS data ready for the processing
load imageData.mat  % imageData matlab
disp(['Start Time:', datestr(now())])
%% Problem Definition

nVar = 4;                                  % Number of Decision Variables
VarMin = -50;                                % Lower Bound of Variables
VarMax = 50;                                % Upper Bound of Variables
MaxIt = 10;                                % Maximum Number of Iterations 1000 25
nPop = 5;                                 % Population Size (Swarm Size) 30 10
runs = 30;                                  % Runs  30 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters 算法
empty_particle.Time= zeros(1, runs);
empty_particle.Methods=[];
empty_particle.BestSol.Position=zeros(runs, nVar);
empty_particle.BestSol.Cost=zeros(1, runs);
empty_particle.BestSol.subPosition=zeros(runs, nVar);
empty_particle.BestSol.subCost=zeros(1, runs);
empty_particle.Curve=zeros(runs,MaxIt);
Methods = {'IAOA'; 'AOA'; 'PSO'; 'ACO'; 'DE'; 'GA'; 'GWO'; 'WOA'; 'EO'; 'AHA'; 'AO'};
Parameters = {...
    @(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction) IAOA(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction);...
    @(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction) AOA(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction); ...
    @(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction) PSO(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction); ...
    @(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction) ACO(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction); ...
    @(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction) DE(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction); ...
    @(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction) GA(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction); ...
    @(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction) GWO(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction); ...
    @(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction) WOA(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction); ...
    @(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction) EO(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction); ...
    @(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction) AHA(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction);...
    @(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction) AO(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction);};
dataNum = imageData.num * numel(imageData.satelliteType) * numel(imageData.Thematic);

startNum = 1;
EndNum = 100;

for i = startNum : EndNum
    t0 = tic;
    currentNum = 0;
    folder  = ['DataProcessor/data/'  num2str(i) ];
    if ~exist(folder,'dir')
        mkdir(folder)
    end
    for satelliteType = 1: numel(imageData.satelliteType)
        for Thematic = 1: numel(imageData.Thematic)
            for num = 1: imageData.num
                model = struct;
                if satelliteType == 1
                    MS = imageData.QuickBird(num, Thematic).MS;
                    PAN = imageData.QuickBird(num, Thematic).PAN;
                else
                    MS = imageData.Gaofen(num, Thematic).MS;
                    PAN = imageData.Gaofen(num, Thematic).PAN;
                end
                [model.MSWV_US, model.F_P, model.F_R, model.F_G, model.F_B, model.F_NIR, ...
                    model.P, model.I, model.bandCoeffs, model.MSWV_DG, model.PANWV_DS, model.W_R, model.W_B, ...
                    model.W_G, model.W_NIR] = ...
                    initInput(MS, PAN);
                model.MS = MS;
                CostFunction=@(x) ERGAS_Index(x, ...
                    model.MSWV_US, model.W_R, model.W_G, model.W_B, ...
                    model.W_NIR, model.F_P, model.F_R, ...
                    model.F_G, model.F_B, model.F_NIR, ...
                    model.P, model.I, model.bandCoeffs, model.MSWV_DG);
                
                data=repmat(empty_particle,numel(Parameters),1);
                for run = 1 : runs
                    for j = 1 : numel(Parameters)
                        t1 = tic;
                        [BestSol,Curve]=Parameters{j,1}(nPop, MaxIt, VarMin, VarMax, nVar, CostFunction);
                        data(j).Methods = Methods{j,1};
                        data(j).BestSol.Position(run, :) = BestSol.Position;
                        data(j).BestSol.Cost(run) = BestSol.Cost;
                        data(j).BestSol.subPosition(run, :) = BestSol.subPosition;
                        data(j).BestSol.subCost(run) = BestSol.subCost;
                        data(j).Curve(run,:) = Curve;
                        data(j).Time(run)=toc(t1);
                    end
                end
                PATH = [folder '/P' num2str(mod(currentNum, dataNum)+1)];
                save(PATH, 'data', 'model');
                currentNum = currentNum + 1;
            end
        end
    end
    % 数据处理
    DataProcessorResult = DataProcessor(imageData, Methods, folder);
    PATH = ['DataProcessor/DataProcessorResult/R' num2str(i)];
    save(PATH, 'DataProcessorResult');
    disp(['Running： ', num2str((i/EndNum)*100), '%, in timing: ', num2str(toc(t0))]);
end
disp(['End Time:', datestr(now())])