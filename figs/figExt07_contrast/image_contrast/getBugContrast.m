function [bugContrast,fishContrast,fileNames]=getBugContrast 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate object contrast, user selects the object
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize variables 
global BIGEYEROOT
filepath=[BIGEYEROOT, '/data/vision/imagecontrast/'];
    close all;	% Close all figure windows except those created by imtool.
    imtool close all;	% Close all figure windows created by imtool.

    fileNames={[filepath,'Giant Centipede.jpeg'],...
        [filepath,'brownCentipede.jpg'],...
        [filepath,'centipedes.jpg'],...
        [filepath,'maxresdefault.jpg'],...
        [filepath,'maxresdefault2.jpg'],...
        [filepath,'maxresdefault3.jpg'],...
        [filepath,'millipede.jpg'],...
        [filepath,'millipedeWood.jpg'],...
        [filepath,'redBlackMillipede.jpg'],...
        [filepath,'ScolonpendraPolymorp.jpg'],...
        [filepath,'blackCentipede.png'],...
        [filepath,'VietnameseCentipede.jpg'],...
        [filepath,'002.jpg'],...
        [filepath,'amazon_cichlid.png'],...
        [filepath, 'amazon_fish.png'],...
        [filepath,'arowanawild.jpg'],...
        [filepath,'fishContrast.png']};
    warning('Do not close figure window during execution!')
    for i=1:length(fileNames)
        baseFileName=fileNames{i};
        fullFileName=fullfile(baseFileName);

        if~exist(fullFileName,'file')
            fullFileName=baseFileName;
            if~exist(fullFileName,'file')
                errorMessage=sprintf('Error: %s does not exist in the search path.',fullFileName);
                uiwait(warndlg(errorMessage));
                return;
            end
        end

        figure();

        I=imread(fullFileName);
        imshow(I,[]);
        axis on;
        title('Select the object of interest by clicking on the figure and dragging the mouse until the object is enclosed');
        set(gcf,'Position',get(0,'Screensize'));

        hFHObject=imfreehand();

        binaryImage=hFHObject.createMask();
        grayImage=rgb2gray(I);
        labeledImage=bwlabel(binaryImage);
        blackMaskedImage=grayImage;
        blackMaskedImage(~binaryImage)=0;
        measurements=regionprops(binaryImage,grayImage,'area','Perimeter');
        area=measurements.Area;
        perimeter=measurements.Perimeter;

        structBoundaries=bwboundaries(binaryImage);
        xy=structBoundaries{1};
        x=xy(:,2);
        y=xy(:,1);

        subplot(2,2,1);
        imshow(grayImage,[]); axis on; title('Original Image')
        drawnow;
        subplot(2,2,2)
        imshow(binaryImage); axis on; title('Binary Masked Image');
        drawnow;
        subplot(2,2,1)
        hold on; plot(x,y,'LineWidth',2); title('Superimposed Perimeter')
        drawnow;

        topLine=min(x); bottomLine=max(x);
        leftColumn=min(y); rightColumn=max(y);
        width=bottomLine-topLine+1; height=rightColumn-leftColumn+1;
        croppedImage=imcrop(blackMaskedImage,[topLine,leftColumn,width,height]);

        subplot(2,2,3);
        imshow(croppedImage); axis on; title('Cropped Image to Selection')
        drawnow;

        LObject=croppedImage(croppedImage~=0);
        LmaxObject=max(LObject); LmeanObject=mean(LObject);

        insideMasked=grayImage;
        insideMasked(binaryImage)=0;

        subplot(2,2,4);
        imshow(insideMasked); axis on; title('Background Image, Removed Selection')
        drawnow;

        LBackground=insideMasked(insideMasked~=0);
        LmaxBackground=max(LBackground); LmeanBackground=mean(LBackground);

        if i<13
            bugContrast(i)=(LmeanObject-LmeanBackground)/LmeanBackground;
        else
            fishContrast(i-12)=(LmeanObject-LmeanBackground)/LmeanBackground;
        end

        close all;
    end

save([BIGEYEROOT 'figExt07_contrast/image_contrast/imageContrastValues'],'bugContrast','fishContrast','fileNames');