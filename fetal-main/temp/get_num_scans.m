function scans = get_num_scans()
E_num = dir('E*');
[~,slist] = system(sprintf('exam_info %s',E_num.name));
slis = split(slist,'images');
num_scans=0;

for i =1:length(slis)
    words= split(slis{i},' ');
    for j=1:length(words)
        word=words{j};
        if strcmp(word,'RTSSFSE')
            %if str2num(words{end-1}) > (num_scans)
                temp_num = str2num(words{end-1});
                num_scans=num_scans +temp_num ;
            %end
        end
    end
end

scans = (num_scans);