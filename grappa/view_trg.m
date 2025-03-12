% visualizing the target points dynamically
tst=zeros(dx,dy,dz);
for i=1:90:size(trg,2)
    tst(trgx(i),trgy(i),trgz(i))=1;
    im(squeeze(tst(4,:,:)))
    pause(0.5);
end