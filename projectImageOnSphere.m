function projectImageOnSphere(filename,radius)
% projectImageOnSphere(filename,radius) projects an image (preferable a
% surface map of for example a planet) on a sphere with a given radius.
% 
% Inputs:
% filename      Filename of the image of the surface map. Call imformats to
%               see the supported image formats.
% radius        The radius of the sphere. This can be used to scale the
%               sphere such that is can be applied for a specific
%               application (e.g., the radius of a planet such that you can
%               plot the trajectory around it)
% 
% 
% 
% Author:           Jaimy Hess
% E-mail:           J.R.Hess@student.tudelft.nl
%                   jaimyhess@gmail.com
% Date modified:    June 5th, 2015
% 
    % Load the file
    c=imread(filename);
  
    % Initalize longitude and colatitude range
    longitude = deg2rad(linspace(0,360,size(c,2)));
    colatitude = deg2rad(linspace(90,-90,size(c,1)));
    
    % Initialize empty variables
    X = zeros(size(c,1),size(c,2));
    Y = zeros(size(c,1),size(c,2));
    Z = zeros(size(c,1),size(c,2));
    % Convert spherical to cartesian coordinates
    for j = 1:length(colatitude)
        for i = 1:length(longitude)
            [X(j,i), Y(j,i), Z(j,i)] = sph2cart( ...
                                            longitude(i),...
                                            colatitude(j),...
                                            radius);
        end
    end
    % Print surface
    surface(X,Y,Z,c,'EdgeColor','none','FaceColor','texturemap');
    axis equal
    
end