function disp=track2disp(tracks)
%d=track2disp(tracks)
% createss a handy structure for bead positions and displacements useful in
% TFM
%tracks: output of track.m [x,y, t, id] or [x,y,z,t,id]
%output:  astructure d, containing all positions and displacements for each
%time.  
%disp(t).r: positions at time t,  
%disp(t).dr: displacements at time t relative to t=1
%NOTE:
% this code rejects all particles that are not present at all time points.
%CREATED BY ERD January 2009
%KNOW ISSUES:
% in some cases, you might want to keep some particles that are present at
% only some timepoints

[r, c] = size(tracks);
if c==4
    d=2;
else
    d=3;
end

ti=min(tracks(:,d+1));
%make first time 1
if ti~=1
    tracks(:,d+1)=tracks(:,d+1)-ti+1;
end
tf=max(tracks(:,d+1));%找出最多追踪的帧数

np = max(tracks(:,d+2));

%remove all tracks that are not present for the whole dataset
ndel=0;
for id=1:np
    ind = find(tracks(:,d+2)==id);
    if length(ind)<(tf);
        tracks(ind,:)=[];
        ndel=ndel+1;%如果第id个粒子追踪帧数少于最大帧数，则将其删除，并记录删除的粒子处ndel
    end
end
display(['deleting ',int2str(ndel),'/',int2str(np) ' incomplete trajectories']);

for t=1:tf
    tks = tracks(tracks(:,d+1)==t,:);
    if t==1
        r0=tks(:,1:d);
    end
   disp(t).r = tks(:,1:d);
   disp(t).dr = tks(:,1:d)-r0;
end
