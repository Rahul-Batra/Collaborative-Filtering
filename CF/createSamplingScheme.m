
function samplingLocations = createSamplingScheme (sizeImage, patternName, patternParameter)

% assume that the center of kspace is at the center of the image

switch patternName
 case 'random'
  
  samplingLocations = rand (sizeImage) < patternParameter;
  samplingLocations = + samplingLocations;
  
 case 'random_axis_aligned_lines'
  
  % num. of full lines = num. of half lines / 2
  patternParameter = patternParameter / 2;
  
  % 33% of lines in the center of kspace (densely sampled)
  % rest of the lines are randomly taken on each "side",
  % but an equal number of lines are taken on each side
  centerLines = round (patternParameter * 0.33);
  if mod (centerLines, 2) == 1
    centerLines = centerLines + 1;
  end
  sparseLines = round (patternParameter - centerLines);
  if mod (sparseLines, 2) == 1
    sparseLines = sparseLines + 1;
  end
  
  %
  
  samplingLocations = zeros(sizeImage);
  
  samplingLocations(sizeImage(1)/2-centerLines/2+1:sizeImage(1)/2+centerLines/2,:) = 1;
  
  T = randperm(sizeImage(1)/2-centerLines/2);
  
  samplingLocations(T(1:sparseLines/2),:) = 1;
  
  T = randperm(sizeImage(1)/2-centerLines/2);
  
  samplingLocations(sizeImage(1)/2+centerLines/2+T(1:sparseLines/2),:) = 1;
  
  % take the center of kspace to the corners
  samplingLocations = fftshift (samplingLocations);
  
 case 'radial_approx_on_cartesian'
  
  numOfLines = patternParameter;
  
  center = (sizeImage + 2) / 2;
  
  halfDiagonal = norm (sizeImage) / 2;
  
  stepLength = 0.5;
  
  numOfSteps = round (halfDiagonal / stepLength + 1);
  
  samplingLocations = zeros (sizeImage);
  
  for lineNum = 0 : numOfLines - 1
    
    theta = 2 * pi * lineNum / numOfLines;
    
    direction = [cos(theta) sin(theta)];
    
    for stepNum = 0 : numOfSteps - 1
      
      location = round (center + direction * stepNum * stepLength);
      
      if ( (location(1) >= 1) && (location(1) <= sizeImage(1)) && (location(2) >= 1) && (location(2) <= sizeImage(2)) )
        
        samplingLocations (location(1), location(2)) = 1;
        
      end
      
    end
    
  end
  
  % take the center of kspace to the corners
  samplingLocations = fftshift (samplingLocations);
  
 case 'radial'
  
  numOfLines = patternParameter;
  
  stepLength = min ([2*pi/sizeImage(1) 2*pi/sizeImage(2)]);
  
  numOfSteps = round (norm ([2*pi 2*pi]) / stepLength + 1);
  
  samplingLocations = [0 0];
  
  count = 2;
  
  for lineNum = 0 : numOfLines - 1
    
    theta = 2 * pi * lineNum / numOfLines;
    
    direction = [cos(theta) sin(theta)];
    
    for stepNum = 1 : numOfSteps - 1
      
      location = direction * stepNum * stepLength;
      
      if norm (location) > pi
        break;
      end
      
      if ( (location(1) >= -pi) && (location(1) < pi) && (location(2) >= -pi) && (location(2) < pi) )
        
        samplingLocations (count, :) = location;
        
        count = count + 1;
        
      else
        
        break;
        
      end
      
    end
    
  end
  
 case 'radial_jitterDirection_uniformLine_approx_on_cartesian'
  
  numOfLines = patternParameter;
  
  center = (sizeImage + 2) / 2;
  
  halfDiagonal = norm (sizeImage) / 2;
  
  stepLength = 0.5;
  
  numOfSteps = round (halfDiagonal / stepLength + 1);
  
  samplingLocations = zeros (sizeImage);
  
  temp = linspace (0, 2 * pi, numOfLines + 1); % temporary variable storing the start and end points of each of the (numOfLines) intervals
  
  for i = 1:length(temp)-1
    fullTheta(i) = temp(i) + rand*(temp(i+1) - temp(i)); % jitter sampled angles stored as vector
  end
  
  for lineNum = 1 : numOfLines
    
    theta = fullTheta(lineNum);
    
    direction = [cos(theta) sin(theta)];
    
    for stepNum = 0 : numOfSteps - 1
      
      location = round (center + direction * stepNum * stepLength);
      
      if ( (location(1) >= 1) && (location(1) <= sizeImage(1)) && (location(2) >= 1) && (location(2) <= sizeImage(2)) )
        
        samplingLocations (location(1), location(2)) = 1;
        
      end
      
    end
    
  end
  
  % take the center of kspace to the corners
%   samplingLocations = fftshift (samplingLocations);
  
 case 'radial_jitterDirection_uniformLine_approx_on_cartesian_symmetric'
  
  if mod (patternParameter, 2) == 1
    numOfLines = patternParameter + 1;
  else
    numOfLines = patternParameter;
  end
  
  center = (sizeImage + 2) / 2;
  
  halfDiagonal = norm (sizeImage) / 2;
  
  stepLength = 0.5;
  
  numOfSteps = round (halfDiagonal / stepLength + 1);
  
  samplingLocations = zeros (sizeImage);
  
  temp = linspace (0, pi, round (numOfLines/2) + 1); % temporary variable storing the start and end points of each of the (numOfLines) intervals
  
  for i = 1:length(temp)-1
    half1Theta(i) = temp(i) + rand*(temp(i+1) - temp(i)); % jitter sampled angles stored as vector
  end
  half2Theta = half1Theta + pi;
  fullTheta = [half1Theta half2Theta];
  
  for lineNum = 1 : numOfLines
    
    theta = fullTheta(lineNum);
    
    direction = [cos(theta) sin(theta)];
    
    for stepNum = 0 : numOfSteps - 1
      
      location = round (center + direction * stepNum * stepLength);
      
      if ( (location(1) >= 1) && (location(1) <= sizeImage(1)) && (location(2) >= 1) && (location(2) <= sizeImage(2)) )
        
        samplingLocations (location(1), location(2)) = 1;
        
      end
      
    end
    
  end
  
  % take the center of kspace to the corners
  samplingLocations = fftshift (samplingLocations);
  
 case 'radial_jitterDirection_jitterLine_approx_on_cartesian'
  
  numOfLines = patternParameter;
  
  center = (sizeImage + 2) / 2;
  
  stepLength = 2;
  
  minRadius = sizeImage(1)/2;
  
  samplingLocations = zeros (sizeImage);
  
  temp = linspace(0,2*pi,numOfLines+1); % temporary variable storing the start and end points of each of the (numOfLines) intervals
  
  for i = 1:length(temp)-1
    fullTheta(i) = temp(i) + rand * (temp(i+1) - temp(i)); % jitter sampled angles stored as vector
  end
  
  for lineNum = 1 : numOfLines
    
    theta = fullTheta (lineNum);
    direction = [cos(theta) sin(theta)];
    
    % find the length of the line in the particular direction
    if (theta == 0) || (theta == pi) || (theta == pi/2) || (theta == 1.5*pi)
      dist1 = minRadius;
    end
    
    if ((0 < theta) && (theta < pi/4)) || ((0.75*pi < theta) && (theta < 1.25*pi)) || ((1.75*pi < theta) && (theta < 2*pi))
      dist1 = floor(minRadius/abs(cos(theta)));
    end
    
    if ((pi/4 < theta) && (theta < 0.75*pi)) || ((1.25*pi < theta) && (theta < 1.75*pi))
      dist1 = floor(minRadius/abs(sin(theta)));
    end
    
    if (theta == 0.25*pi) || (theta == 0.75*pi) || (theta == 1.25*pi) || (theta == 1.75*pi)
      dist1 = floor(sqrt(2)*minRadius);
    end
    
    totalLineInterval = round(dist1/stepLength); % Given the step length and the angle, find the total number of intervals
    
    for i = 1:totalLineInterval
      
      line(i) = stepLength * (i-1) + rand * stepLength;
      
      location = round (center + direction * line(i));
      
      if ( (location(1) >= 1) && (location(1) <= sizeImage(1)) && (location(2) >= 1) && (location(2) <= sizeImage(2)) )
        samplingLocations (location(1), location(2)) = 1;
      end
    end
    
  end
  
  % take the center of kspace to the corners
  samplingLocations = fftshift (samplingLocations);
  
 case 'radial_jitterDirection_jitterLine_approx_on_cartesian_symmetric'
  
  if mod(patternParameter,2)
    numOfLines = patternParameter+1;
  else
    numOfLines = patternParameter;
  end
  
  center = (sizeImage + 2) / 2;
  
  stepLength = 2;
  
  minRadius = sizeImage(1)/2;
  
  samplingLocations = zeros (sizeImage);
  
  temp = linspace(0,pi,round(numOfLines/2)+1); % temporary variable storing the start and end points of each of the (numOfLines) intervals
  
  for i = 1:length(temp)-1
    half1Theta(i) = temp(i) + rand * (temp(i+1) - temp(i)); % jitter sampled angles stored as vector
  end
  half2Theta = half1Theta + pi;
  fullTheta = [half1Theta half2Theta];
  
  for lineNum = 1 : numOfLines
    theta = fullTheta(lineNum);
    direction = [cos(theta) sin(theta)];
    
    % find the length of the line in the particular direction
    if (theta == 0) || (theta == pi) || (theta == pi/2) || (theta == 1.5*pi)
      dist1 = minRadius;
    end
    
    if ((0 < theta) && (theta < pi/4)) || ((0.75*pi < theta) && (theta < 1.25*pi)) || ((1.75*pi < theta) && (theta < 2*pi))
      dist1 = floor(minRadius/abs(cos(theta)));
    end
    
    if ((pi/4 < theta) && (theta < 0.75*pi)) || ((1.25*pi < theta) && (theta < 1.75*pi))
      dist1 = floor(minRadius/abs(sin(theta)));
    end
    
    if (theta == 0.25*pi) || (theta == 0.75*pi) || (theta == 1.25*pi) || (theta == 1.75*pi)
      dist1 = floor(sqrt(2)*minRadius);
    end
    
    totalLineInterval = round(dist1/stepLength); % Given the step length and the angle, find the total number of intervals
    
    for i = 1:totalLineInterval
      
      line(i) = stepLength * (i-1) + rand * stepLength;
      
      location = round (center + direction * line(i));
      
      if ( (location(1) >= 1) && (location(1) <= sizeImage(1)) && (location(2) >= 1) && (location(2) <= sizeImage(2)) )
        samplingLocations (location(1), location(2)) = 1;
      end
    end
    
  end
  
  % take the center of kspace to the corners
  samplingLocations = fftshift (samplingLocations);
  
 case 'radial_uniformDirection_jitterLine_approx_on_cartesian'
  
  numOfLines = patternParameter;
  
  center = (sizeImage + 2) / 2;
  
  stepLength = 2;
  
  minRadius = sizeImage(1)/2;
  
  samplingLocations = zeros (sizeImage);
  
  fullTheta = 2*pi/numOfLines:2*pi/numOfLines:2*pi; % Uniform division of angles stored in a vector
  
  for lineNum = 1 : numOfLines
    
    theta = fullTheta(lineNum);
    direction = [cos(theta) sin(theta)];
    
    % find the length of the line in the particular direction
    if (theta == 0) || (theta == pi) || (theta == pi/2) || (theta == 1.5*pi)
      dist1 = minRadius;
    end
    
    if ((0 < theta) && (theta < pi/4)) || ((0.75*pi < theta) && (theta < 1.25*pi)) || ((1.75*pi < theta) && (theta < 2*pi))
      dist1 = floor(minRadius/abs(cos(theta)));
    end
    
    if ((pi/4 < theta) && (theta < 0.75*pi)) || ((1.25*pi < theta) && (theta < 1.75*pi))
      dist1 = floor(minRadius/abs(sin(theta)));
    end
    
    if (theta == 0.25*pi) || (theta == 0.75*pi) || (theta == 1.25*pi) || (theta == 1.75*pi)
      dist1 = floor(sqrt(2)*minRadius);
    end
    
    totalLineInterval = round(dist1/stepLength); % Given the step length and the angle, find the total number of intervals
    
    for i = 1:totalLineInterval
      
      line(i) = stepLength*(i-1) + rand*stepLength;
      
      location = round (center + direction * line(i));
      
      if ( (location(1) >= 1) && (location(1) <= sizeImage(1)) && (location(2) >= 1) && (location(2) <= sizeImage(2)) )
        samplingLocations (location(1), location(2)) = 1;
      end
    end
    
  end
  
  % take the center of kspace to the corners
  samplingLocations = fftshift (samplingLocations);
  
 case 'radial_uniformDirection_jitterLine_approx_on_cartesian_symmetric'
  
  numOfLines = patternParameter;
  
  center = (sizeImage + 2) / 2;
  
  stepLength = 2.0;
  
  minRadius = sizeImage(1)/2;
  
  samplingLocations = zeros (sizeImage);
  
  half1Theta = 2*pi/numOfLines:2*pi/numOfLines:pi; % Uniform division of angles stored in a vector
  half2Theta = half1Theta + pi;
  fullTheta = [half1Theta half2Theta];
  
  for lineNum = 1 : numOfLines
    
    theta = fullTheta(lineNum);
    direction = [cos(theta) sin(theta)];
    
    % find the length of the line in the particular direction
    if (theta == 0) || (theta == pi) || (theta == pi/2) || (theta == 1.5*pi)
      dist1 = minRadius;
    end
    
    if ((0 < theta) && (theta < pi/4)) || ((0.75*pi < theta) && (theta < 1.25*pi)) || ((1.75*pi < theta) && (theta < 2*pi))
      dist1 = floor(minRadius/abs(cos(theta)));
    end
    
    if ((pi/4 < theta) && (theta < 0.75*pi)) || ((1.25*pi < theta) && (theta < 1.75*pi))
      dist1 = floor(minRadius/abs(sin(theta)));
    end
    
    if (theta == 0.25*pi) || (theta == 0.75*pi) || (theta == 1.25*pi) || (theta == 1.75*pi)
      dist1 = floor(sqrt(2)*minRadius);
    end
    
    totalLineInterval = round(dist1/stepLength); % Given the step length and the angle, find the total number of intervals
    
    for i = 1:totalLineInterval
      
      line(i) = stepLength*(i-1) + rand*stepLength;
      
      location = round (center + direction * line(i));
      
      if ( (location(1) >= 1) && (location(1) <= sizeImage(1)) && (location(2) >= 1) && (location(2) <= sizeImage(2)) )
        samplingLocations (location(1), location(2)) = 1;
      end
    end
    
  end
  
  % take the center of kspace to the corners
  samplingLocations = fftshift (samplingLocations);
  
 otherwise
  
  error ('error in patternParameter in createSamplingScheme');
  
end

return
