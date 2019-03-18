%% Automated detection and quantification of Babesiosis infection
% Babesiosis is a malaria-like, tick-borne parasitic disease common in
% parts of the US, and present sporadically throughout the rest of the
% world. Our goal here is to EXPLORE options for developing an automated
% function to DETECT the presence of the Babesia parasite in thin blood
% smears, and to QUANTIFY the portion of RBCs in a sample that are
% infected.
%
% Brett Shoelson,PhD; Avi Nehemiah

%% Credits
% All images used in this demo can be found
% <http://www.cdc.gov/dpdx/babesiosis/gallery.html here>. We gratefully
% acknowledge the CDC's Division of Parasitic Diseases and Malaria (DPDM).
% Images are used by permission.
%

%% OUR FOCUS
%
% Exploration
% UI usage
% Segmentation
% Pre-processing

%% Clean slate
clear; close all;clc; 

%% Babesiosis Images
babesiosisDir = fullfile(pwd,'.\BloodSmearImages\babesiosis');
% Convenient referencing of a collection of images? 
imgSet = imageSet(babesiosisDir) %#ok<*NOPTS>
methods(imgSet)

%% Create a display of all Babesiosis images
togglefig('Babesiosis Images')
ax = gobjects(imgSet.Count,1);
for ii = 1:imgSet.Count
	ax(ii) =...
		subplot(floor(sqrt(imgSet.Count)),ceil(sqrt(imgSet.Count)),ii);
	[~,currName] = fileparts(imgSet.ImageLocation{ii});
	imshow(read(imgSet,ii))
	title([num2str(ii),') ' currName],...
		'interpreter','none','fontsize',7)
end
expandAxes(ax);

%% To develop an algorithm, let's select a "Target Image"
targetImgNum = 4;
togglefig('Babesiosis Images')
[~,imName] = fileparts(imgSet.ImageLocation{targetImgNum});
set(ax,'xcolor','r','ycolor','r',...
	'xtick',[],'ytick',[],'linewidth',2,'visible','off')
set(ax(targetImgNum),'visible','on');

%% 
% targetImage = imread(imgSet.ImageLocation{targetImgNum});
% OR:
% targetImage = read(imgSet,targetImgNum);
% OR:
targetImage = getimage(ax(targetImgNum));
togglefig('Target Image')
clf
imshow(targetImage)
title(imName,'interpreter','none','fontsize',12);

%% Segmentation
% When we "segment" an image, we distinguish the regions of interest (ROIs)
% from the non-ROI portion, generally creating a binary mask of what we
% want to qualify, quantify, track, etc. Segmentation is a critical part of
% many image processing problems, and is worth considering in some depth.
% Let's select an image for development of a segmentation algorithm:
%
% NOTE:
%   First, we have to decide _what_ we want to segment...The cells? The
%   infections? Let's initially target the cells, and then determine which
%   of them are infected.
%
%   There are many approaches...
%   Let's explore using an image segmentation app
% 
% colorThresholder(targetImage)
% imageSegmenter(targetImage) %-> 'segmentImageFcn.m'

%% Apply auto-generated segmenter

cellMask = segmentImageFcn(targetImage);
togglefig('Cell Mask')
imshow(cellMask);

%% Improving the result
% doc('Remove small objects image')
cellMask = bwareaopen(cellMask,100);
togglefig('Cell Mask',true)
imshow(cellMask);

%% Looks like a good start! But now we need to separate CONTIGUOUS cells
% Segmenting touching objects is one of the most challenging tasks in image
% processing.
%  Let's consider a few approaches:
%  * Using edges
%  * Using watershed
%  * Using other information? Shape?

%% Try detecting edges
% doc edge
edges = edge(rgb2gray(targetImage));
togglefig('Edge Mask')
subplot(1,2,1)
imshow(targetImage)
subplot(1,2,2)
imshow(edges);

%% 
% segmentImage(rgb2gray(targetImage)) %Using Brett's app
edges = edge(rgb2gray(targetImage),'LOG',0.001);
togglefig('Edge Mask',1)
imshow(edges);

%% Improving the result
edges = bwareaopen(edges,60,8);
togglefig('Edge Mask')
imshow(edges);

%% Combine the edges (logically) with the segmented regions
togglefig('Cell Mask')
tmp = cellMask & ~edges;
tmp = bwareaopen(tmp,100);
imshow(tmp);

%% Improve the edge mask?
% Unfortunately, we haven't yet separated the contiguous cells!
%
% imageMorphology(edges)

morphed1 = imclose(edges, strel('Disk',3,4));
morphed1 = bwmorph(morphed1, 'skeleton', Inf);
morphed1 = bwmorph(morphed1, 'spur', Inf);
togglefig('Edge Mask',true)
imshow(morphed1);

%% Consider further "cleaning" of the mask: 
% imageRegionAnalyzer(morphed1);
morphed1 = bwpropfilt(morphed1,'Perimeter',[80 Inf]);
togglefig('Edge Mask',true)
imshow(morphed1);

%
imshow(cellMask & ~morphed1);

%% Beautiful!
cellMask = cellMask & ~imdilate(morphed1,strel('disk',1));
togglefig('Target Image',true)%Automatically clear the figure
tmpAx(1) = axes('units','normalized',...
	'position',[0.025 0.525 0.25 0.25]);
imshow(targetImage);
tmpAx(2) = axes('units','normalized',...
	'position',[0.025 0.225 0.25 0.25]);
imshow(label2rgb(bwlabel(cellMask), @jet, 'k', 'shuffle'))
tmpAx(3) = axes('units','normalized',...
	'position',[0.3 0.025 0.675 0.95]);
imshow(targetImage)
showMaskAsOverlay(0.4,cellMask,'g')
expandAxes(tmpAx);

%% Science vs Art?
% ...But the truth of the matter is that the more we customize on one
% image, the less general our approaches may become. Let's automate the
% analysis of the images using this algorithm:
%
% (Create 'refinedMask.m')

%% 
togglefig('Babesiosis Images',true)
refreshImages
for ii = 1:imgSet.Count
	mask = refinedMask(getimage(ax(ii)));
	showMaskAsOverlay(0.5,mask,'b',[],ax(ii))
	drawnow
end
expandAxes(ax);

%% ALTERNATIVE APPROACHES TO CONSIDER
% So where does that leave us?
% 
% The fact of the matter is, when objects are contiguous, getting a good
% segmentation can be difficult. But obtaining a good segmentation mask can
% be essential to getting good analysis results.
%
% Here are a couple of alternatives to consider:

%% Watershed Segmentation
% doc('segment touching regions images')
%
togglefig('Exploration',true) %Automatically clear the figure
grayscale = rgb2gray(targetImage);
imshow(grayscale)
%
wsImg = watershed(grayscale);
showMaskAsOverlay(0.8, wsImg==0, 'g')
title('Classic "Oversegmentation"','fontsize',14);

%% How can we suppress those small catchment basins?

doc('suppress minima image')

%% Improving that result:
% Consider the image as a "topography"
togglefig('Cell Mask',true)
clf
ax3(1) = subplot(1,2,1);
ds = 3;
surf(im2double(grayscale(1:ds:end,1:ds:end)));
shading interp;
rotate3d on
set(gca,'view',[0 90],...
	'xlim',[0 100],'ylim',[0 100],'zlim',[0.3 1])
title('Original Grayscale','fontsize',12)
drawnow
% IMHMIN:
newGrayscale = imhmin(grayscale,13);
ax3(2) = subplot(1,2,2);
surf(im2double(newGrayscale(1:ds:end,1:ds:end)));
shading interp
rotate3d on
set(gca,'view',[0 90],...
	'xlim',[0 100],'ylim',[0 100],'zlim',[0.3 1])
title('Grayscale with Suppressed Minima','fontsize',12)
linkprop(ax3,'view');
colormap(flipud(parula));

%%
togglefig('Exploration')
newGrayscale = imhmin(grayscale,13);
wsImg = watershed(newGrayscale);
showMaskAsOverlay(1,wsImg == 0, 'g');
title('Better...','fontsize',14);

%% So now we can use this calculation to improve our segmentation mask
togglefig('Cell Mask',true)
cellMask = segmentImageFcn(targetImage);
cellMask = bwareaopen(cellMask,30);
imshow(cellMask);

%% Impose watershed lines
% Just like we did with the edges calculated earlier, we can break our
% regions along watershed lines
togglefig('Cell Mask',true)
wsEdges = wsImg == 0;
wsEdges = bwareaopen(wsEdges,200,8);
cellMask(wsEdges) = 0;
imshow(cellMask);

%% But again...this doesn't necessarily generalize well:
togglefig('Babesiosis Images')
for ii = 1:imgSet.Count
	tmpMask = refinedMask2(imgSet,ii);
	showMaskAsOverlay(0.5,tmpMask,'g',true,ax(ii))
	drawnow
end
expandAxes(ax);

% Further reading on watershed segmentation:
% http://www.mathworks.com/company/newsImgletters/articles/the-watershed-transform-strategies-for-image-segmentation.html
% http://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/

%% Using the "shapes" of the objects of interest...
% Consider treating the RBCs as circles
% doc imfindcircles
% circleFinder(grayscale)

%%
detectCircles = @(x) imfindcircles(x,[20 35], ...
	'Sensitivity',0.89, ...
	'EdgeThreshold',0.04, ...
	'Method','TwoStage', ...
	'ObjectPolarity','Dark');
[centers, radii, metric] = detectCircles(grayscale);

togglefig('Target Image',true)
imshow(targetImage)
viscircles(centers,radii,'edgecolor','b')
title(sprintf('%i Cells Detected',numel(radii)),'fontsize',14);

%% Again, we check to see how robust the approach is:
togglefig('Babesiosis Images')
for ii = 1:imgSet.Count
	[centers,radii] = detectCircles(rgb2gray(read(imgSet,ii)));
	delete(findall(ax(ii),'type','line'))
	viscircles(ax(ii),centers,radii,'edgecolor','b')
	drawnow
end
expandAxes(ax);

%% We're almost there! Now we need only quantify infection...
% We want to differentiate cells suspected of being infected from otherwise
% healthy cells. Masking and evaluating the regions might be useful. But
% where is the useful information?

imtool(grayscale);

%% However...the variability in the images is still troubling
imtool close all
%
togglefig('Exploration',1)
infectionThreshold = 135;
infection = grayscale <= infectionThreshold;
subplot(1,2,1)
imshow(targetImage)
showMaskAsOverlay(1,infection,'r');
subplot(1,2,2)
tmpImg = read(imgSet,16);
imshow(tmpImg)
infection = rgb2gray(tmpImg) <= infectionThreshold;
showMaskAsOverlay(1,infection,'r');

%% One last step
% doc imhistmatch
tempImage = read(imgSet,16);
matchedImage = imhistmatch(tempImage,targetImage);
togglefig('Exploration',true)
subplot(2,1,1)
imshow(targetImage);title('target')
subplot(2,1,2)
imshowpair(tempImage,matchedImage,'montage');
title('Image 16;                     Image 16, HistMatched');

%% So now...
imtool(rgb2gray(imfuse(tempImage,matchedImage,'montage')));

%% So which cells, and what fraction of cells, are infected?
imtool close all
%
[centers,radii] = detectCircles(grayscale);
isInfected = false(numel(radii),1);
nCells = numel(isInfected);

%
% Creating a "mesh" can be useful:
x = 1:size(grayscale,2);
y = 1:size(grayscale,1);
[xx,yy] = meshgrid(x,y);
% xx(1:5,1:5),yy(1:5,1:5)
togglefig('Grayscale',true);
imshow(grayscale)
infectionMask = false(size(grayscale));
for ii = 1:numel(radii)
	mask = hypot(xx - centers(ii,1), yy - centers(ii,2)) <= radii(ii);
	currentCellImage = grayscale;
	currentCellImage(~mask) = 0;
	infection = ...
		currentCellImage > 0 & currentCellImage < infectionThreshold;
	% OPTIONALLY set a range for the calculated regionproperties of
	%   infections:
	% infection = bwpropfilt(infection,'Eccentricity',[0,0.925]);
	infectionMask = infectionMask | infection;
	isInfected(ii) = any(infection(:));
	if isInfected(ii)
		showMaskAsOverlay(0.3,mask,'g',[],false);
	end
end
showMaskAsOverlay(0.5,infectionMask,'r',[],false)
title(sprintf('%i of %i (%0.1f%%) Infected',...
	sum(isInfected),numel(isInfected),...
	100*sum(isInfected)/numel(isInfected)),...
	'fontsize',14,'color','r');

%% Is this more generalizable?
togglefig('Babesiosis Images',true)
refreshImages;
drawnow
%
for ii = 1:imgSet.Count
	[pctInfection,centers,radii,isInfected,infectionMask] = ...
		testForInfection(getimage(ax(ii)),targetImage,...
		infectionThreshold,detectCircles);
	title(ax(ii),...
		['Pct Infection: ', num2str(pctInfection,2),...
		' (' num2str(sum(isInfected)),...
		' of ' num2str(numel(isInfected)) ')']);
	viscircles(ax(ii),centers,radii,'edgecolor','b')
	% createCirclesMask
	infectedCellsMask = createCirclesMask(targetImage,...
		centers(isInfected,:),...
		radii(isInfected));
	showMaskAsOverlay(0.3,infectedCellsMask,'g',ax(ii),false);
	showMaskAsOverlay(0.5,infectionMask,'r',ax(ii),false);
	drawnow
end
expandAxes(ax);

%% Now generate a report
% options = {'html','doc','latex','ppt','xml','pdf'};
% myDoc = publish('ParasitologyDemo.m',options{6})
% winopen(myDoc)
return %For publishing purposes...

%% And then, a note on speeding up analyses using parallel computing
%gcp; %Get/start a "parpool"
tic; %#ok
pctInfection = zeros(imgSet.Count,1);
for ii = 1:imgSet.Count
	img = imread(imgSet.ImageLocation{ii});
	[pctInfection(ii),centers,radii,isInfected,infectionMask] = ...
		testForInfection(img,targetImage,...
		infectionThreshold,detectCircles);
end
tSerial = toc

tic
pctInfection = zeros(imgSet.Count,1);
allImages = [imgSet.ImageLocation];
parfor ii = 1:imgSet.Count
	img = imread(allImages{ii});
	[pctInfection(ii),centers,radii,isInfected,infectionMask] = ...
		testForInfection(img,targetImage,...
		infectionThreshold,detectCircles);
end
tParallel = toc

%% < MACHINE LEARNING/CLASSIFICATION >
%
% BACK TO POWERPOINT...
%
clear;close all;clc; %#ok

%% Can we differentiate types of infections from blood-smear images?
% Create an image set representing Motorbikes and airplanes. These are
imgSet = imageSet(fullfile(pwd,'.\BloodSmearImages'),...
	'recursive')  
disp(['Your imageSet contains ', num2str(sum([imgSet.Count])),...
	' images from ' num2str(numel(imgSet)) ' classes.']);

% imageSetViewer(imgSet)

%% What do these images look like?
subset = select(imgSet,1:3);
subsetNames = [subset.ImageLocation];
subsetLabels = {};
for ii = 1:numel(subset)
	subsetLabels{ii} = repelem({subset(ii).Description},subset(ii).Count,1);%#ok
end
subsetLabels = vertcat(subsetLabels{:});
togglefig('Sample Images',1)
[hpos,hdim] = distributeObjects(numel(subset),0.05,0.95,0.01);
[vpos,vdim] = distributeObjects(3,0.95,0.05,0.025);
ax = gobjects(numel(subset),1);
[hind,vind] = meshgrid(1:numel(imgSet),1:subset(1).Count);
for ii = 1:numel(subsetNames)
	ax(ii) = axes('Units','Normalized',...
		'Position',...
		[hpos(hind(ii)) vpos(vind(ii)) hdim vdim]);
	imshow(subsetNames{ii});
	title(subsetLabels{ii},'fontsize',8)
end
expandAxes(ax);

%% First we PARTITION the imageSet into training and test sets
[trainingSets, testSets] = partition(imgSet,0.7,'randomized');

%% The we create a visual BAG OF FEATURES to describe the training set:
bag = bagOfFeatures(trainingSets);
                                                                                          editorwindow;
%% Visualize Feature Vectors 
togglefig('Encoding',true)
for ii = 1:numel(imgSet)
	img = read(imgSet(ii), randi(imgSet(ii).Count));
	featureVector = encode(bag, img);
	subplot(numel(imgSet),2,ii*2-1);
	imshow(img);
	title(imgSet(ii).Description)
	subplot(numel(imgSet),2,ii*2);
	bar(featureVector);
	set(gca,'xlim',[0 bag.VocabularySize])
	title('Visual Word Occurrences');
	if ii == numel(imgSet)
		xlabel('Visual Word Index');
	end
	if ii == floor(numel(imgSet)/2)
		ylabel('Frequency of occurrence');
	end
end
                                                                                                      editorwindow;
%% TRAIN category classifier on the training set
classifier = trainImageCategoryClassifier(trainingSets,bag);

%% EVALUATE the classifier on the test-set images:
[confMat,knownLabelIdx,predictedLabelIdx,predictionScore] = ...
	evaluate(classifier,testSets);
avgAccuracy = mean(diag(confMat));
togglefig('Prediction')
imagesc(confMat)
colorbar

%% Now we can use the classifier to PREDICT class membership!
togglefig('Prediction')
ii = randi(size(imgSet,2));
img = read(imgSet(ii),randi(imgSet(ii).Count));
[labelIdx, predictionScore] = predict(classifier,img);
bestGuess = classifier.Labels(labelIdx);
actual = imgSet(ii).Description;
imshow(img)
t = title(['Best Guess: ',bestGuess{1},'; Actual: ',actual]);                                                          editorwindow
if strcmp(bestGuess{1},actual)
	set(t,'color',[0 0.7 0])
else
	set(t,'color','r')
end

%% We can easily try other classifiers using the classificationLearner
% Here we recreate the bagOfFeatures from all images, and cast it to a
% table to facilitate working with the classificationLearner app
bag = bagOfFeatures(imgSet)
infectionData = double(encode(bag, imgSet));
InfectionImageData = array2table(infectionData);
infectionType = categorical(repelem({imgSet.Description}',...
	[imgSet.Count], 1));
InfectionImageData.infectionType = infectionType;                                                                           editorwindow;

%% Use the new features to train a model and assess its performance using 
classificationLearner

%% NOT AS GOOD AS YOU'D LIKE? TRY AGAIN WITH A 2-CLASS MODEL & A CUSTOM EXTRACTOR:
clear('trainedClassifier*')
% First we create a single-source of non-plasmodium images:
nonPlasmodium = imgSet(~strcmp({imgSet.Description},'plasmodium'));
nonPlasmodium = imageSet([nonPlasmodium.ImageLocation]);
nonPlasmodium.Description = 'nonPlasmodium';

% And then a Plasmodium-nonPlasmodium one:
twoClassSet = cat(1,nonPlasmodium,...
	imgSet(strcmp({imgSet.Description},'plasmodium')))

% Re-Learn Bag of Feature representation
extractorFcn = @customParasitologyFcn;
% We will use:
%  * a CUSTOM EXTRACTOR,
%  * ALL FEATURES
%  * LARGER VOCABULARY

bag = bagOfFeatures(twoClassSet,...
	'CustomExtractor',extractorFcn,...
	'StrongestFeatures',1,...
	'VocabularySize',1000);

infectionData = double(encode(bag, twoClassSet));
InfectionImageData = array2table(infectionData);
infectionType = categorical(repelem({twoClassSet.Description}',...
	[twoClassSet.Count], 1));
InfectionImageData.infectionType = infectionType;
classificationLearner

% EXPORT MODEL -> trainedClassifier
%%
trainedClassifier

%% Now we can use the classifier to PREDICT class membership!
% PREDICT the infection type in a randomly selected test image
ii = randi(size(twoClassSet,1));
jj = randi(twoClassSet(ii).Count);
img = read(twoClassSet(ii),jj);
%img = read(imgSet(ii),randi(imgSet(ii).Count));
togglefig('Test Image'); set(gcf,'color','w');
imshow(img)
% Add code here to invoke the trained classifier
imagefeatures = double(encode(bag, img));
% Find two closest matches for each feature
[bestGuess, predictionScore] = predict(trainedClassifier,imagefeatures);
% Display the string label for img
if strcmp(char(bestGuess),twoClassSet(ii).Description)
	titleColor = [0 0.8 0];
else
	titleColor = 'r';
end
title(sprintf('Best Guess: %s;\nActual: %s',...
	char(bestGuess),twoClassSet(ii).Description),...
	'color',titleColor)
editorwindow;

%% STILL not as good as you'd like?
% How can we improve the predictive value?
% 
% * Better images
% * More images!
% * Preprocessing
% * Non-gridded point selection
% * Custom extractor (see: bagOfFeatures)
% * More images
% * More images
% * More images
% * More images

%% Links to auxiliary tools referenced herein
%
% <http://www.mathworks.com/matlabcentral/fileexchange/47956 editorwindow>
%
% <http://www.mathworks.com/matlabcentral/fileexchange/18220 togglefig>
%
% <http://www.mathworks.com/matlabcentral/fileexchange/18291 expandAxes>
%
% <http://www.mathworks.com/matlabcentral/fileexchange/19706 ExploreRGB>
%
% <http://www.mathworks.com/matlabcentral/fileexchange/22108 showMaskAsOverlay>
%
% <http://www.mathworks.com/matlabcentral/fileexchange/48859 Image Segmenter app>
%
% <http://www.mathworks.com/matlabcentral/fileexchange/23697 Image Morphology app>
%
% <http://www.mathworks.com/matlabcentral/fileexchange/34365 Circle Finder app>
%
% <http://www.mathworks.com/matlabcentral/fileexchange/47905 createCirclesMask>