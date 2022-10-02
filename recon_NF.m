% This scripts runs reconstruction from speckle data on our examples. 
% Note: the code is runing on a GPU

% members in data_param:
%       is_PSI: using phase shifting interferometry or not. If yes, group
%       every three captured image to get interfered term.
%       image_num_use: number of images used for reconstruction
%       angle_m_list: angle and slope used for compensating phase ramp
%
% members in alg_param:
%       fw: width of blurring kernel, to remove the DC from the input intensity images.
%       T_tau: the size of local window beeing corelated. By defualt this should be 
%               equivalent to the speckle spread, but due to 
%               computational considerations other sizes can be used. 
%       T_delta: the support around a window were correlation is computed.
%           Ideally T_delta/2 correspond to maximal displacment where we expect some ME
%           correlation. Can be set lower to speed computation.     
%       max_sft: the desired size of the output image is max_sft x max_sft (must be squared image)
%       iter: number of iteration to update the algorithms


function recon_img = recon_NF(frame_all,data_param,alg_param)


mean_spk_sprd = alg_param.T_tau;
if length(mean_spk_sprd)==1
    tg=[1:ceil(mean_spk_sprd)]; tg=tg-mean(tg(:));
    [gx,gy]=meshgrid(tg,tg); 
    gr=sqrt(gx.^2+gy.^2);
    mean_spk_sprd=(gr<=floor(mean_spk_sprd/2));
end

% interval between local correlation is coputed
stp_size=floor(alg_param.T_tau/2);

[n1,n2]=size(frame_all{1,1});


%define 2 windows: Wwin0, the size of the local window we correlate in the
%speckle image. Wwin, the local window in the target image. 
%Following eq(24) in the [11], 
%the window from the target image should be wider,
%proportional to the speckle spread from one source.
Wwin0=ones(alg_param.T_tau,alg_param.T_tau);
Wwin=conv_fft2(Wwin0,mean_spk_sprd);
if alg_param.norm_window 
    Wwin = Wwin/max(Wwin(:));
end

Wwin0=padarray(Wwin0,(size(Wwin)-size(Wwin0))/2);
bdr_size=(min(n1,n2)-alg_param.max_sft)/2;

%phase shifting interferometery funciton as in eq(18) int he paper
PSI_fun = @(I) (double(I{1}) + double(I{2})*exp(2j*pi/3) + double(I{3})*exp(4j*pi/3)); 

cu_all = {};

if data_param.is_PSI
    p_id_list = 1:data_param.image_num_use/3;        
    l = size(frame_all,2)/3;
else
    p_id_list = 1:data_param.image_num_use;
    l = size(frame_all,2);    
end


I = frame_all{1,1};
[wy,wx] = size(I);
[X,Y] = meshgrid((1:wx)-floor((wx+1)/2), (1:wy)-floor((wy+1)/2));


for p_id = p_id_list
    p_id
    image_all = {};    
    for i = 1:size(frame_all,1)   
        if data_param.is_PSI           
            img = PSI_fun(frame_all(i,[p_id,p_id+l,p_id+l*2]));                                      
        else
            img = double(frame_all{i,p_id });
        end
        img = img/65536;
        
        image_all = [image_all, img]; 
    end
    image_sum = sum(cat(3,image_all{:}),3);

    % remove DC component
    if alg_param.fw ~= -1
        if(isreal(image_sum))
            blurred_image_sum = imgaussfilt(image_sum,alg_param.fw);
        else
            blurred_image_sum = imgaussfilt(real(image_sum),alg_param.fw)+1j*imgaussfilt(imag(image_sum),alg_param.fw);
        end
        image_sum = image_sum-blurred_image_sum;
    end
    
    I =  image_sum;
    
    % compensate phase ramp
    if isfield(data_param,'angle_m_list')
        rot_X = imrotate(X,data_param.angle_m_list(1,p_id),'nearest','crop');
        %rot_X = imrotate(X,-data_param.angle_m_list(1,p_id),'nearest','crop');
        ramp_correction = exp(1j*rot_X*data_param.angle_m_list(2,p_id));
    else
        ramp_correction = ones(size(I));
    end
    
    % calculate local correlations
    tic
    cu=getWinCorrW(I,stp_size,alg_param.T_delta,Wwin0,ramp_correction );
    toc
        
    % merge local correlations from different observations
    if p_id == 1
        cu_sum = cu;
    else
        cu_sum = cu_sum + cu;
    end
    
    cu_sum =  real(cu_sum);  
end
cu_sum = cu_sum.*(cu_sum>0);        
cu_sum = cu_sum/length(p_id_list);


cu = cu_sum;

sc=0.1;
maxItr2=1; maxItr1=alg_param.iter;

disp('begin optimizeWinCorW')
datetime('now')	

% reconstruction target
targ = sc*gpuArray(cu/max(cu(:)));

%optimization
[res_local,m_err]=optimizeWinCorW(size(I)-bdr_size*2,stp_size,alg_param.T_delta,gpuArray(Wwin),targ,bdr_size,maxItr1,maxItr2,alg_param.suppress_DC);

recon_img = gather(real(res_local));

end
