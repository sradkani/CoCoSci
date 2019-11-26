function motifs = generateMotifs(alphabet, maxMotifLength)    
% preallocate motif cell
    motifs = cell(1, maxMotifLength);

    % generate motifs
    for i = 1:maxMotifLength
        % each entry in motifs will be a matrix containing all motifs of a
        % given length
        motifs{i} = permn(alphabet, i);

        % remove motifs that can be described by shorter motifs

        % this is the case when longer motifs repeat themselves
        % e.g. ACAC can be represented by 'AC' motif alone

        % find partitions of motifs of length i 
        % (e.g. i = 4 -> mod(4, 0:3) -> 4 0 0 1 -> 
        % find 0's to get valid partitions)
        partitions = find(mod(i, 1:i-1) == 0); 

        % iterate through partition, e.g. if i = 4, partitions are 1 and 2
        for partition = partitions
            % duplicate 1:partition columns of matrix to check 
            % which rows contain duplications
            duplicatedMat = repmat(motifs{i}(:,1:partition), [1 i/partition]);

            % rows that are in duplicated matrix are redundant
            [~, redundantRows] = intersect(motifs{i}, duplicatedMat, 'rows');

            % delete redundant rows
            motifs{i}(redundantRows, :) = [];

        end
    end