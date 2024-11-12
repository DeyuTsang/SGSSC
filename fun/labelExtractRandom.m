function out = labelExtractRandom(gnd, labelNum, labelRatio)


if ~exist('gnd','var')
    error('LabelError',' Input hypergraph must have labels...');
end

vnum = length(gnd);

if (isempty(labelNum))
    if (~isempty(labelRatio))
        labelNum = ceil(vnum*labelRatio);
    else
        error('LabelError',' at least one of labelNum and labelRatio must be non-empty...');
    end
end

if (labelNum >= vnum)
    error('LabelError','HAT:label number must be smaller than the number of instances...');
end

% initial assignment [1,...,1,0,...,0]
out = sparse(ones(1,labelNum), 1:labelNum, true, 1, vnum);

% randomly permute out until all labels appear in out
out = out(randperm(vnum));   

count = 0;
while (~checkLabels(gnd, out))
    out = out(randperm(vnum));
    count = count + 1;
    if (count >= 1000)
        printLabels(hypergraph);
        error('HAT:labelError','HAT:impossible to make all labels appear in labelMask, try to increase labelNum. Maybe some label has very few vertices. Remove those labels if it is the case (use compactHypergraphKeepLabel)...');
    end
end

% reinitialize the result (because there might be a bug with matlab's sparse matrix...)
out = full(out);
out = sparse(out);


function checkResult = checkLabels(gnd, out)
        % check if all labels appear in out
        checkResult = true;
        for i=1:length(gnd)
            if (isempty(find(out & gnd(i,:),1)))
                checkResult = false;
                break;
            end
        end



