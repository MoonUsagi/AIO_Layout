mImg = imread('Input.jpg');
mImg = rgb2gray(mImg);
mImg = im2bw(mImg);
figure;
imshow(mImg);
