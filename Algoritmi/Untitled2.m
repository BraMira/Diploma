k = 4;
l = 3;
h = 6;

for i = 1:2,
    for j = 1:1
        if k>l
            k = k+l
            if k>h
                k = k+h
                break
            else
                disp('no')
                continue
            end
        else
            disp('nonon')
            continue
        end
    end
end

disp(k)

for i = 1:2
    if k>0
        disp('done')
    end
end
            
            