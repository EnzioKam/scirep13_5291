foldername = 'data';
mkdir(foldername);

a = 10;
dims = [10, 30, 50, 100];
cno_list = [4];
reps = 100;

for cno = cno_list
    for d = dims
        if cno == 1 || cno == 4 || cno == 20
            tmatrix1 = zeros(d, reps);
        elseif cno == 2
            tmatrix1 = zeros(d, reps);
            rmatrix1 = zeros(d, d, reps);
        elseif cno == 5
            tmatrix1 = zeros(d, reps);
            rmatrix1 = zeros(d, d, reps);
            rmatrix2 = zeros(d, d, reps);
        end
        
        for i = 1:reps
            if cno == 1 || cno == 4 || cno == 20
                tmatrix1(:, i) = make_translate(a, d);
            elseif cno == 2
                tmatrix1(:, i) = make_translate(a, d);
                rmatrix1(:, :, i) = make_rotate(d);
            elseif cno == 5
                tmatrix1(:, i) = make_translate(a, d);
                rmatrix1(:, :, i) = make_rotate(d);
                rmatrix2(:, :, i) = make_rotate(d);
            end
        end
        
        filename = [foldername '/C' num2str(cno) '_' num2str(d)];
        if cno == 1 || cno == 4 || cno == 20
            save(filename, 'a', 'd', 'cno', 'reps', 'tmatrix1')
        elseif cno == 2
            save(filename, 'a', 'd', 'cno', 'reps', 'tmatrix1', 'rmatrix1')
        elseif cno == 5
            save(filename, 'a', 'd', 'cno', 'reps', 'tmatrix1', 'rmatrix1', 'rmatrix2')
        end
    end
end