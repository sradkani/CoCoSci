function curve = randomXCurve(alphabet, maxMotifLength, sequence, maxRandomX, delta, alpha)
% takes as input the sequence and parameters for the inferring HMM and 
% returns a curve of random(X) for each point in the  sequence 
% (starting from the 2nd state)

% this is a vector that will contain randomness at each point
curve = nan(1, length(sequence)-1);

% compute random(X) at each point in the sequence
for i = 1:length(sequence)-1
    curve(i) = findRandomness(alphabet, maxMotifLength, sequence(1:i+1), delta, alpha);
    
    curve(i) = curve(i); %%./maxRandomX(i);

end


    