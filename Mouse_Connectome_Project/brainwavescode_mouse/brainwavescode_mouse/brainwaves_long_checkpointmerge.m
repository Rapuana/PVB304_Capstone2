function brainwaves_long_checkpointmerge(indir)
% merge brainwaves_long checkpoint files

if ~strcmp(indir(end),filesep)
    indir=[indir filesep];
end

list=dir([indir '*END*']);
list={list.name};

for k=1:length(list)
    basename=list{k}(1:end-7);
    
    END=load([indir basename 'END.mat']);
    ncheckpoints=END.checkpoint;
    
    in=load([indir basename '1.mat']);
    
    
    time=in.time(1):in.outdt:END.time(end);
    
    soln=zeros(size(in.soln,1),length(time));
    
    soln(:,1:size(in.soln,2))=in.soln;
    filledto=size(in.soln,2);
    
    overlap=in.ictime/in.outdt + 1;
    
    for j=2:ncheckpoints-1
        in=load([indir basename num2str(j) '.mat']);
        newfilledto=filledto+size(in.soln,2)-overlap;
        soln(:,filledto+1:newfilledto)=in.soln(:,overlap+1:end);
        filledto=newfilledto;
    end
    
    newfilledto=filledto+size(END.soln,2)-overlap;
    soln(:,filledto+1:newfilledto)=END.soln(:,overlap+1:end);
    filledto=newfilledto;
    
    if filledto~=length(time)
        error('soln and time are different lengths')
    end
    
    out=END;
    out.soln=soln;
    out.time=time;
    
    save([indir basename(1:end-11)],'-struct','out')
    
end


% delete
for k=1:length(list)
    basename=list{k}(1:end-7);
    delete([indir basename '*.mat'])
end