gpuDevice(1);

filepath = './';
filelist = {'main_cam_150um','main_cam_150um2','main_cam_150um3','main_cam_150um4','main_cam_150um5',...
    'main_cam_150um6','main_cam_150um7'};

%filepath = './data/speckles1108/';
%filelist = {'virus_175','neuron_175_2','virus_125_2','neuron_125'};


%parpool(3);
%recon_imgs = {};

% suppressDC = false;
% 
% for f_iter = 1:2
%     load([filepath,filelist{f_iter}]);    
%     image_num_use = 48;
%     for sample_rate = [1,2,4]
%         is_PSI = true;
%         if is_PSI
%             p_id_list = 1:image_num_use/3;
%             l = size(frame_all,2)/3;
%         else
%             p_id_list = 1:image_num_use;
%             l = size(frame_all,2);    
%         end
%         
%         recon_img = recon_NF([filepath,filelist{f_iter}], is_PSI,p_id_list,l,sample_rate,suppressDC);
%         save(sprintf('recon/local_recon_1017/%s_%dimgs_%dsr.mat',filelist{f_iter},image_num_use,sample_rate),'recon_img');            
%     end
% end
suppress_DC = 2;
more_phase = false;

sample_rate = 1;
for f_iter = 7
    load([filepath,filelist{f_iter}],'frame_all');    
    if f_iter ==9
        single_recon=true;        
        frame_all_temp = reshape(frame_all,[4,54]);
        frame_all = cell(1,54);
        for i=1:54
            frame_all{1,i} = sum(cat(3,frame_all_temp{:,i}),3);
        end   
    elseif f_iter ==10
        single_recon=true;                
    else       
        single_recon=true;        
        frame_all = frame_all(:,55);
        
%         single_recon=false;        
%         frame_all = frame_all(:,1:54);
%         frame_all = reshape(frame_all,[3,18])';
%         frame_all = reshape(frame_all,[1,54]);       
    end
    
    if single_recon
        image_num_use = 1;%54;
        frame_all = frame_all(:,1:image_num_use);                 
        is_PSI = false;          
        fw = 7;   
        more_phase = false;
    else
        more_phase = false;      
        if more_phase      
            image_num_use = 108;       
        else
            image_num_use = 54;                   
        end                   
        frame_all = frame_all(:,1:image_num_use);         
        
        is_PSI = true;
        fw = 37;
    end
    
    %m_base = 0.016;
    %image_num_use = 1;
    
    alpha = -0.005;
    m_base = alpha*(6.5/20)* (2*pi/0.58)* 1.625; % alpha * length_per_pixel* k* |dt|;
    
    for m_base = 0.016
        w_size = 25;
        data_params = struct('is_PSI',is_PSI,'image_num_use', image_num_use);
        %alg_params = struct('fw',fw,'w_size',49,'max_delta',201,'max_sft',201,'iter',500,'suppress_DC',suppress_DC,'norm_window',false);
        alg_params = struct('fw',fw,'w_size',w_size,'max_delta',151,'max_sft',401,'iter',200,'suppress_DC',suppress_DC,'norm_window',true,'more_phase',more_phase);
              
        
        if ~single_recon
            data_params.angle_m_list = [[0:30:150,0:15:165]; [m_base*ones(1,6),2*m_base*ones(1,12)]];
        end
        
        recon_img = recon_NF(frame_all, data_params, alg_params);
        save_dir = sprintf('recon/local_recon_0311/%s',filelist{f_iter});
        
        if ~exist( save_dir,'dir')
            mkdir(save_dir);
        end
        
        save_filename = sprintf('%s/img%d_sr%d_fw%d_ws%d_md%d_ms%d_iter%d',save_dir, image_num_use,sample_rate,...
        alg_params.fw, alg_params.w_size, alg_params.max_delta, alg_params.max_sft, alg_params.iter);
        
        if suppress_DC == 1
            save_filename = [save_filename,'_SD'];
        elseif suppress_DC == 2
            save_filename = [save_filename,'_SD2'];            
        end
        
        if data_params.w_neg
            save_filename = [save_filename,'_N'];
        end        
        
        if alg_params.norm_window
            save_filename = [save_filename,'_nw'];            
        end
        
        if alg_params.more_phase
            save_filename = [save_filename,'_mp'];            
        end        

        if isfield(data_params,'angle_m_list')
            save_filename = [save_filename,sprintf('_mbase_%.3f', m_base);];                        
        end
                
        save([save_filename,'.mat'],'recon_img');            
        save_img([save_filename,'.png'],recon_img);
    end
end
