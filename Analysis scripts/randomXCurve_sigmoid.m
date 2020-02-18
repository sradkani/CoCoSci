function curve = randomXCurve_sigmoid(alphabet, maxMotifLength, sequence, delta, alpha)
% takes as input the sequence and parameters for the inferring HMM and 
% returns a curve of random(X) for each point in the  sequence 
% (starting from the 2nd state)

% this is a vector that will contain randomness at each point
curve = nan(1, length(sequence)-1);

% compute random(X) at each point in the sequence
for i = 1:length(sequence)-1
    
    randomness = findRandomness(alphabet, maxMotifLength, sequence(1:i+1), delta, alpha);
    
    curve(i) = 1/(1+exp(-0.2*randomness));
    
end


    