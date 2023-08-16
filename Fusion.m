
function FusionResult = Fusion(Sol, MSWV_US, W_R, W_G, W_B, W_NIR, F_P, F_R, F_G, F_B, F_NIR, P, I, bandCoeffs, MS_ORG)

row = size(Sol,1);
empty_particle.Result=[];

FusionResult=repmat(empty_particle,row,1);
Sol=(Sol+50)/100;
for i = 1 : row
    FusionResult(i).Result(:,:,1) = MSWV_US(:,:,1) + W_R.*(Sol(i,1)*F_P + (1-Sol(i,1))*F_R).*(P-I);
    FusionResult(i).Result(:,:,2) = MSWV_US(:,:,2) + W_G.*(Sol(i,2)*F_P + (1-Sol(i,2))*F_G).*(P-I);
    FusionResult(i).Result(:,:,3) = MSWV_US(:,:,3) + W_B.*(Sol(i,3)*F_P + (1-Sol(i,3))*F_B).*(P-I);
    FusionResult(i).Result(:,:,4) = MSWV_US(:,:,4) + W_NIR.*(Sol(i,4)*F_P + (1-Sol(i,4))*F_NIR).*(P-I);
    
    for j = 1:size(MS_ORG, 3)
        FusionResult(i).Result(:,:,j) = FusionResult(i).Result(:,:,j)*bandCoeffs(j);
    end
end


end