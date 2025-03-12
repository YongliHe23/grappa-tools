% visualizing the source points dynamically
tst=zeros(dx,dy,dz);
for i=1:90:size(src,2)
    tst(srcx(i),srcy(i),srcz(i))=1;
    im(squeeze(tst(3,:,:)))
    pause(0.5);
end