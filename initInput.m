function [MSWV_US, F_P, F_R, F_G, F_B, F_NIR, P, I, bandCoeffs, MSWV_DG, PANWV_DS, W_R, W_B, W_G, W_NIR] = initInput(MS, PAN)

%% Make the PAN and MS data ready for the processing ʹPAN��MS�������ô���׼��
PAN = PAN(:,:, 1);
MSWV_db  = double(MS);
PANWV_db = double(PAN);
% MS_ORG   = double(imageData{1}.MS);

%% Resizing, Upsampling the MS data to the size of PAN ������С����MS�����ϲ�����PAN�Ĵ�С

MSWV_US  = imresize(MSWV_db,  1/4, 'bicubic');
MSWV_US  = imresize(MSWV_US,  4,   'bicubic');
% MSWV_US  = MSWV_db;
MSWV_DG  = MSWV_US;
% MSWV_US = MSWV_db;
% MSWV_DG  = MSWV_US;
if size(PANWV_db, 1) == size(MSWV_US, 1)
    PANWV_DS = PANWV_db;
end
if size(PANWV_db, 1) ~= size(MSWV_US, 1)
    PANWV_DS = imresize(PANWV_db,  1/4, 'bicubic');
end
% PANWV_DS = PANWV_db;
% PANWV_DS = imresize(PANWV_db,  1/4, 'bicubic');

%% Seperating the spectral bands �����״�

R   = MSWV_US(:,:,1); 
G   = MSWV_US(:,:,2); 
B   = MSWV_US(:,:,3); 
NIR = MSWV_US(:,:,4); 

%% Data Normialization ���ݹ淶��
num = size(MSWV_US,3);
bandCoeffs = zeros(1, num);
for i=1:num
    bandCoeffs(i)      = max(max(MSWV_US(:,:,i)));
    MSWV_US(:,:,i)     = MSWV_US(:,:,i)/bandCoeffs(i);
end

P = PANWV_DS;
panCoeff = max(max(max(P)));
P = P/panCoeff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pansharpening Framework ���񻯿��

% Edge detecotrs 

lamda = 10^-9;
eps   = 10^-10;

% PAN weight in edge detector

F_P   = expEdge(P, lamda,eps);

% MS weight in edge detector
Red_W = max(max(R));      % Red component
R     = R/Red_W;
F_R   = expEdge(R, lamda, eps);

Green_W = max(max(G));    % Green component
G       = G/Green_W;
F_G     = expEdge(G, lamda,eps);

Blue_W  = max(max(B));    % Blue component
B       = B/Blue_W;
F_B     = expEdge(B, lamda,eps);

NIR_W=max(max(NIR));      % NIR component
NIR=NIR/NIR_W;
F_NIR = expEdge(NIR, lamda,eps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extracting primitive detail map ��ȡ������ϸ��ͼ

W = impGradDes(MSWV_US,P);   % Optimal weights for spectral bands

I     = W(1)*MSWV_US(:,:,1) + W(2)*MSWV_US(:,:,2) + W(3)*MSWV_US(:,:,3) + W(4)*MSWV_US(:,:,4); 
P     = (P-mean(P(:)))*std(I(:))/std(P(:)) + mean(I(:));  % Histogram matching

W_R   = 4*R./(R + B + G + NIR);
W_B   = 4*B./(R + B + G + NIR);
W_G   = 4*G./(R + B + G + NIR);
W_NIR = 4*NIR./(R + B + G + NIR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end