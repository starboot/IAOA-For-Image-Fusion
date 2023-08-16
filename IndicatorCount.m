
function tmpN = IndicatorCount(ObjectiveAssesmentData, IAOAr)
tmpN = zeros(1, IAOAr);

for kk = 1 : IAOAr
    temp = ObjectiveAssesmentData(kk, :);
    nn = 0;
    for zz = 1 : IAOAr
        if kk ~= zz
            if temp.ERGAS < ObjectiveAssesmentData(zz, :).ERGAS
                nn = nn + 1;
            end
            if temp.SAM < ObjectiveAssesmentData(zz, :).SAM
                nn = nn + 1;
            end
            if temp.RASE < ObjectiveAssesmentData(zz, :).RASE
                nn = nn + 1;
            end
            if temp.RMSE < ObjectiveAssesmentData(zz, :).RMSE
                nn = nn + 1;
            end
            if temp.UIQI > ObjectiveAssesmentData(zz, :).UIQI
                nn = nn + 1;
            end
            if temp.CC > ObjectiveAssesmentData(zz, :).CC
                nn = nn + 1;
            end
        end
    end
    tmpN(kk) = nn;
end
end