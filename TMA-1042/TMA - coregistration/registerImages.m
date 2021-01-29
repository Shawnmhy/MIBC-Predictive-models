function [MOVINGREG] = registerImages(MOVING,FIXED)
%registerImages  Register grayscale images using auto-generated code from Registration Estimator app.
%  [MOVINGREG] = registerImages(MOVING,FIXED) Register grayscale images
%  MOVING and FIXED using auto-generated code from the Registration
%  Estimator app. The values for all registration parameters were set
%  interactively in the app and result in the registered image stored in the
%  structure array MOVINGREG.

% Auto-generated by registrationEstimator app on 25-May-2020
%-----------------------------------------------------------


% Convert RGB images to grayscale
FIXED = rgb2gray(FIXED);
MOVING = rgb2gray(MOVING);

% Default spatial referencing objects
fixedRefObj = imref2d(size(FIXED));
movingRefObj = imref2d(size(MOVING));

% Intensity-based registration
[optimizer, metric] = imregconfig('multimodal');
metric.NumberOfSpatialSamples = 500;
metric.NumberOfHistogramBins = 50;
metric.UseAllPixels = true;
optimizer.GrowthFactor = 1.050000;
optimizer.Epsilon = 1.50000e-06;
optimizer.InitialRadius = 6.25000e-03;
optimizer.MaximumIterations = 100;

% Align centers
fixedCenterXWorld = mean(fixedRefObj.XWorldLimits);
fixedCenterYWorld = mean(fixedRefObj.YWorldLimits);
movingCenterXWorld = mean(movingRefObj.XWorldLimits);
movingCenterYWorld = mean(movingRefObj.YWorldLimits);
translationX = fixedCenterXWorld - movingCenterXWorld;
translationY = fixedCenterYWorld - movingCenterYWorld;

% Coarse alignment
initTform = affine2d();
initTform.T(3,1:2) = [translationX, translationY];

% Apply Gaussian blur
fixedInit = imgaussfilt(FIXED,1.000000);
movingInit = imgaussfilt(MOVING,1.000000);

% Normalize images
movingInit = mat2gray(movingInit);
fixedInit = mat2gray(fixedInit);

% Apply transformation
tform = imregtform(movingInit,movingRefObj,fixedInit,fixedRefObj,'rigid',optimizer,metric,'PyramidLevels',3,'InitialTransformation',initTform);
MOVINGREG.Transformation = tform;
MOVINGREG.RegisteredImage = imwarp(MOVING, movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);

% Store spatial referencing object
MOVINGREG.SpatialRefObj = fixedRefObj;

end

