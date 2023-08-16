function DataProcessorResult = DataProcessor(imageData, Methods, folder)



%% Fusion result 融合结果
dataNum = imageData.num * numel(imageData.satelliteType) * numel(imageData.Thematic);
currentNum = 0;
empty_particle.Table=[];
empty_particle.tmpN=0;
empty_particle.data=[];
empty_particle.model=[];
empty_particle.FusionResult=[];
T = repmat(empty_particle,dataNum,1);
for Thematic = 1: numel(imageData.Thematic)
    
    for satelliteType = 1: numel(imageData.satelliteType)
        
        for imageNum = 1: imageData.num
            currentNum = currentNum + 1;
            PATH = [folder '/P' num2str(currentNum)];
            load (PATH);
            MS_ORG = model.MSWV_DG;
            MS = model.MS;
            runs = size(data(1).BestSol.subPosition, 1);
            IAOAr = runs * 2;
            num = numel(data) + IAOAr - 1;
            FusionData = zeros(num, 4);
            for i = 1 : numel(data)
                if i ==  1
                    FusionData(i : runs, :) = data(i).BestSol.Position;
                    FusionData(runs + 1 : IAOAr, :) = data(i).BestSol.subPosition;
                else
                    FusionData(i + IAOAr - 1, :) = data(i).BestSol.subPosition(floor(rand()*runs)+1, :);
                end
                
            end
            FusionResult = Fusion(FusionData, model.MSWV_US, model.W_R, model.W_G, model.W_B, ...
                model.W_NIR, model.F_P, model.F_R, model.F_G, model.F_B, model.F_NIR, ...
                model.P, model.I, model.bandCoeffs, MS_ORG);
            
            [~, MinIndex] = min(data(1).BestSol.Cost);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Objective Assesment of Fusion Results 融合结果的客观评估
            
            % 越大越好: CC, UIQI
            % 越小越好: ERGAS, RASE, SAM, RMSE
            
            ObjectiveAssesmentData = ObjectiveAssesment(double(MS),FusionResult);
            
            tmpN = IndicatorCount(ObjectiveAssesmentData, IAOAr);
            
            [~, indexB]=max(tmpN);
            
            ObjectiveAssesmentData(IAOAr, :) = ObjectiveAssesmentData(indexB, :);
            
            ObjectiveAssesmentData = ObjectiveAssesmentData(IAOAr : end, :);
            
            tmpN = IndicatorCount(ObjectiveAssesmentData, numel(Methods));
            
            T(currentNum).Table = struct2table(ObjectiveAssesmentData, 'RowNames', Methods);
            
            T(currentNum).tmpN = tmpN(1);
            
            T(currentNum).data = data;
            
            T(currentNum).model = model;
            
            T(currentNum).FusionResult = FusionResult(MinIndex).Result(:,:,1:3);
            
        end
    end
    
end
DataProcessorResult = T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save 保存

% save
end