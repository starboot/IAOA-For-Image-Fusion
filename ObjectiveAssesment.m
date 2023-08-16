function ObjectiveAssesmentData = ObjectiveAssesment(I1,I2)
row = numel(I2);

empty_particle.ERGAS=[];
empty_particle.SAM=[];
empty_particle.RASE=[];
empty_particle.RMSE=[];
empty_particle.UIQI=[];
empty_particle.CC=[];

ObjectiveAssesmentData=repmat(empty_particle,row,1);

for i = 1 : row
    ObjectiveAssesmentData(i).ERGAS   = ERGAS(I1, I2(i).Result, 4);
    ObjectiveAssesmentData(i).SAM     = SAM(I1, I2(i).Result);
    ObjectiveAssesmentData(i).RASE    = RASE(I1, I2(i).Result);
    ObjectiveAssesmentData(i).RMSE    = RMSE(I1, I2(i).Result);
    ObjectiveAssesmentData(i).UIQI    = uqi(I1, I2(i).Result);
    ObjectiveAssesmentData(i).CC      = CC(I1, I2(i).Result);
end

end