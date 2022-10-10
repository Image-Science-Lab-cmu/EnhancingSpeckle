% These reconstruction algorithm are adpated from M. Alterman's work 
% Imaging with local speckle intensity correlations: theory and practice
% https://webee.technion.ac.il/people/anat.levin/papers/AltermanTOG2021codeData.tar.gz 

gpuDevice(1);

filepath = './data/';

% three data correspoinding to:
% 1. first row of Fig 7
% 2. first two rows of Fig S7
% 3. Fig S14

filelist = {'spoke','fluorescent_beads','fluorescent_beads_parafilm'};
f_iter = 3;

data_dir = sprintf('./data/%s',filelist{f_iter});

load(data_dir,'frame_all');    
    
% If single_recon= true, it corresponds to M. Alterman's work 
single_recon=false;        
     
if single_recon
    % does not apply LSI, but using single image to do reconstruction only
    image_num_use = 1;              
    frame_all = frame_all(:,55);
    is_LSI = false;
    
    % Gaussian kernel size to reduce DC components
    % as in (4), we need to subtract a low-pass image to remove DC
    % component
    fw = 7;   
else
    % apply LSI, capture 54 images with 18 ramp 
    % and 3 interfered image per ramp    
    image_num_use = 54;                   
    frame_all = frame_all(:,1:54);
    is_LSI = true;
    
    % gaussian kernel size to reduce DC components
    % while theoretically we do not need to reduce DC in LSI
    % empirically apply a smooth filter will improve the results
    if f_iter == 1
        fw = 101;
    else
        fw = 37;
    end
end

data_params = struct(...
    'is_LSI',is_LSI,...    % using LSI or not
    'image_num_use', image_num_use... % number of image used for reconstruction
);

if ~single_recon    
    alpha = -0.003;
    m_base = -alpha*(6.5/20)* (2*pi/0.58)* 1.625; % -alpha * length_per_pixel* k* |dt|; 
    %negative sign: after 4f sign, the image is flipped
    data_params.angle_m_list = [
        [0:30:150,0:15:165]; % direction of the ramp
        [m_base*ones(1,6),2*m_base*ones(1,12)] % distance of the ramp
    ];
end

alg_params = struct(...
    'fw',fw,...             % gaussian kernel size to reduce DC components 
    'T_tau',25,...          % the size of local window beeing corelated
    'T_delta',151,...       % the support around a window were correlation is computed.
    'max_sft',401,...       % the desired size of the output image
    'iter',200,...          % number of iteration to update the algorithms
    'suppress_DC',true...   % supporess DC of correlation during recon or not
);

recon_img = recon_NF(frame_all, data_params, alg_params);
save_dir = sprintf('recon/%s',filelist{f_iter});

if ~exist( save_dir,'dir')
    mkdir(save_dir);
end

save_filename = sprintf('%s/img%d_fw%d_ws%d_md%d_ms%d_iter%d',save_dir, image_num_use,...
alg_params.fw, alg_params.T_tau, alg_params.T_delta, alg_params.max_sft, alg_params.iter);

if alg_params.suppress_DC
    save_filename = [save_filename,'_SD'];
end

if isfield(data_params,'angle_m_list')
    save_filename = [save_filename,sprintf('_mbase_%.3f', m_base);];                        
end

save([save_filename,'.mat'],'recon_img');            
save_img([save_filename,'.png'],recon_img);
