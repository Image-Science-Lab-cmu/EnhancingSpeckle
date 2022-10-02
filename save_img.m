
function save_img(filename,img)
    %img_temp = zeros(200);
    %img_temp(1:size(img,1), 1:size(img,2)) = rot90(img,2);
    %img = img_temp;
    
    set(gcf,'Position', [100 100 200 200]);    
    set(gca,'LooseInset',get(gca,'TightInset'));
        
    imshow(img,[]); 
    colormap(hot);

    axis off;
    set(gca,'units','pixels'); % set the axes units to pixels
    x = get(gca,'position'); % get the position of the axes
    set(gcf,'units','pixels'); % set the figure units to pixels
    y = get(gcf,'position'); % get the figure position
    set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
    set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
    saveas(gcf,filename);        
end


