%Optimize a match between a local correlation of a speckle image and the
%latent one
function [l_map,m_err]=optimizeWinCorW(targ_size,win_stp_size,w_size,Wwin,targ,bdr_size,maxItr1,maxItr2,suppressDC, init_l_map)

stp_size = 1e-2;
m_err    = 1e17;

for itr=1:maxItr2  
    if ~exist('init_l_map','var')
        i_l_map=rand(targ_size)/1000;
        if  strcmp(class(targ),'gpuArray')
            i_l_map=gpuArray(i_l_map);
        end
    else
        i_l_map=init_l_map;
    end

    [o_l_map,err] = ADAMProj(@(x) eval_fg_getWinCorrW_parallelized(x,win_stp_size,w_size,Wwin,targ,bdr_size,suppressDC),i_l_map,stp_size,maxItr1);
    figure, imagesc(gather(real(o_l_map))); colormap(hot); drawnow;
   
    if err<m_err
        m_err=err;
        l_map=o_l_map;
    end
end

