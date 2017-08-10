function figure_hubs(Parcel_params,thresholds,degree,wd,pc)

figure('Position',[1 1 800 1000]);
subplot(1,3,1);
imagesc(degree(Parcel_params.sorti,:),[0 50]);
hline_new(Parcel_params.transitions,'w',2)
set(gca,'YTick',[]);
set(gca,'YTick',Parcel_params.centers,'YTickLabel',Parcel_params.networks);
set(gca,'XTick',thresholds);
title('degree');
xlabel('thresholds');
subplot(1,3,2);
imagesc(wd(Parcel_params.sorti,:),[-1.5 1.5]);
set(gca,'YTick',[]);
xlabel('thresholds');
set(gca,'XTick',thresholds);
hline_new(Parcel_params.transitions,'w',2);
title('within mod degree');
subplot(1,3,3);
imagesc(pc(Parcel_params.sorti,:),[0 0.75]);
set(gca,'XTick',[],'YTick',[]);
hline_new(Parcel_params.transitions,'w',2)
xlabel('thresholds');
set(gca,'YTick',[]);
set(gca,'XTick',thresholds);
title('participation coefficient');

end