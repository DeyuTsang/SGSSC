function show(fea)
[m,n] = size(fea);
faceW = 32; 
faceH = 32; 
numPerLine =26; 
ShowLine = 20; 

Y = zeros(faceH*ShowLine,faceW*numPerLine); 
for i=0:ShowLine-1 
  	for j=0:numPerLine-1 
    	Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW) = reshape(fea(i*numPerLine+j+1,:),[faceH,faceW]); 
  	end 
end 

imshow(Y,[]);colormap(gray);

end