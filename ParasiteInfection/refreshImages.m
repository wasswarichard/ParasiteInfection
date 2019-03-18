ax = gobjects(imgSet.Count,1);
for ii = 1:imgSet.Count
	ax(ii) = subplot(floor(sqrt(imgSet.Count)),ceil(sqrt(imgSet.Count)),ii);
	currName = imgSet.ImageLocation{ii};
	imshow(imread(currName))
	[~,currName] = fileparts(currName);%Strip out extensions
	title([num2str(ii),') ' currName],...
		'interpreter','none','fontsize',7);
end
expandAxes(ax)
