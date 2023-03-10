path = '/Users/massimilianoberardi/Desktop/Spheroids_results/Images/Optical images_sp growth_time_26/12 days/HT1376/';
files = dir(path)
%%
img_anal_res = zeros(10,4);
for i = 3:11
    spheroid = imread(strcat(path,files(i).name));
    spheroidbin = imbinarize(spheroid(:,:,2),0.2896902);%,'adaptive','ForegroundPolarity','dark');

    [~,treshold] = edge(spheroidbin,"sobel");
    fudge = 0.19;
    BWs = edge(spheroidbin,"sobel",treshold*fudge);
    se90 = strel('line',2,90);
    se0 = strel('line',2,0);
    BWfill = imfill(imdilate(BWs,[se90,se0]),"holes");
    BWsingle = bwareaopen(BWfill,9000);
    %smooth a little
    windowSize = 30 ;
    kernel = ones(windowSize) / windowSize ^ 2;
    BWsingle = imbinarize(conv2(single(BWsingle), kernel, 'same'));
    BWsingle = bwareaopen(BWsingle,12000);
    BWsingle = imfill(BWsingle,"holes");
    BWsingle = imerode(BWsingle,[se90,se0]);
    imshowpair(spheroid,BWsingle,'montage')
    %calculations
    perim = regionprops(bwperim(BWsingle),'Perimeter');
    perim = perim(1).Perimeter;
    img_anal_res(i-2,2) = perim;
    img_anal_res(i-2,1) = perim/2/pi;
    img_anal_res(i-2,3) = bwarea(BWsingle);
    img_anal_res(i-2,4) = 4*pi*bwarea(BWsingle)/perim^2;
    pause(0.5)
end