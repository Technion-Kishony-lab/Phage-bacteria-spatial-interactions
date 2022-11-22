%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Display an example of the two-color plaque assay on agar plates %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Printing Supplmentary Figure 5: 2-color plaque assay\n');
f = filesep;
isoNum = {'31','33'};
imgDir = ['..' f 'input_files' f 'twoColorPlqAssay'];
figFolder = 'Figs';
gamma = 1.2;

%% load images and crop
imgs = zeros(3456,5184,3,numel(isoNum));
for iso = 1:numel(isoNum)
    imgR = imread([imgDir f isoNum{iso} '_RFP.CR2']);
    imgY = imread([imgDir f isoNum{iso} '_YFP.CR2']);
    imgs(:,:,:,iso) = cat(3,imgR(:,:,1),imgY(:,:,2),imgR(:,:,1));

end

load([imgDir f 'crop_data.mat'],'outcropX','outcropY','incropX','incropY');

%% example plaques
plqs = [];
plqs(1).p(1,:) = [540 862];
plqs(1).p(2,:) = [643, 1159];
plqs(2).p(1,:) = [1425, 2885];
plqs(2).p(2,:) = [733, 2444];
orient = [1 -1];
%% plot
g = figure(104); clf;
set(gcf,'units','centimeters','position',[1 1 18 10])
imsz = 8.5;
margin_x = 0.2;
margin_y = 0.5;
axes_pos = [margin_x margin_y imsz imsz];
space_x = 0.5;
r = 1000;
arrSz = 150;
plqNum = 1;
for im = 1:numel(isoNum)
    imgR_crop = imgaussfilt(uint8(imgs(outcropY(1):outcropY(2),outcropX(1):outcropX(2),1,im)),20);
    imgY_crop = imgaussfilt(uint8(imgs(outcropY(1):outcropY(2),outcropX(1):outcropX(2),2,im)),20);
    maxR = double(median(maxk(reshape(imgR_crop,[],1),maxMinPxl)))/255;
    maxY = double(median(maxk(reshape(imgR_crop,[],1),maxMinPxl)))/255;
    minR = double(median(mink(reshape(imgY_crop,[],1),maxMinPxl)))/255;
    minY = double(median(mink(reshape(imgY_crop,[],1),maxMinPxl)))/255;
    img_crop = uint8(imgs(outcropY(1):outcropY(2),outcropX(1):outcropX(2),:,im));
    crpAdjImg = imadjust(img_crop,[minR,minY,minR; maxR,maxY,maxR]);

    ax = axes('units','centimeters','position',axes_pos+(im-1)*[space_x+imsz 0 0 0]);
    imshow(crpAdjImg);

    hold on
    rectangle('position',[size(imgR_crop,2)/2-r size(imgR_crop,2)/2-r 2*r 2*r],'curvature',[1 1],'LineStyle','--');
    text(size(imgR_crop,2)/2, size(imgR_crop,1)/2-300, 'Wildtype T7 dilutions','HorizontalAlignment','center');
    text(size(imgR_crop,2)/2, size(imgR_crop,1)/2-1300, 'Evolved sample dilutions','HorizontalAlignment','center');    
    for pl = 1:2
        annotation('textarrow',ax.Position(1)/g.Position(3)+[plqs(im).p(pl,1)-arrSz plqs(im).p(pl,1)]/size(imgR_crop,2)*ax.Position(3)./g.Position(3),...
            ax.Position(2)/g.Position(4)+(size(imgR_crop,1)-[plqs(im).p(pl,2)-orient(pl)*arrSz plqs(im).p(pl,2)])/size(imgR_crop,1).*ax.Position(4)./g.Position(4),...
            'string',num2str(plqNum),'Color','black');
        plqNum = plqNum+1;
    end
end
print([figFolder f 'SuppFigure5'],'-dpng','-r300');
    
    
    
    