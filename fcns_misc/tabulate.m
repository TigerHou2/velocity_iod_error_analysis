%% tabulate.m
function tabulate(dat,format)

if ~exist('format','var')
    format = 1;
end

if matches(format,'%')
    disp('Using percent notation...')
    regexprep(latex(vpa(sym(dat(:,:,1)*100),10)),...
    '([+-]?\d+\.?\d*)','\n${char(compose(''%.2f'',str2num($1)))}\\%')
    regexprep(latex(vpa(sym(dat(:,:,2)*100),10)),...
    '([+-]?\d+\.?\d*)','\n${char(compose(''%.2f'',str2num($1)))}\\%')
    regexprep(latex(vpa(sym(dat(:,:,3)*100),10)),...
    '([+-]?\d+\.?\d*)','\n${char(compose(''%.2f'',str2num($1)))}\\%')
else
    disp('Using scientific notation...')
    regexprep(...
    regexprep(...
    regexprep(latex(vpa(sym(dat(:,:,1)),10)),...
        '([0-9]+\.[0-9]+)','${num2str(str2num($1),''%.3e'')}'),...
        '(e[+-][0-9]+)',' \\times 10^{${strtok($1,''e'')}}'),...
        '((?<=\{[+-])[0-9])','${regexprep($1,''^0*'','''')}')
    regexprep(...
    regexprep(...
    regexprep(latex(vpa(sym(dat(:,:,2)),10)),...
        '([0-9]+\.[0-9]+)','${num2str(str2num($1),''%.3e'')}'),...
        '(e[+-][0-9]+)',' \\times 10^{${strtok($1,''e'')}}'),...
        '((?<=\{[+-])[0-9])','${regexprep($1,''^0*'','''')}')
    regexprep(...
    regexprep(...
    regexprep(latex(vpa(sym(dat(:,:,3)),10)),...
        '([0-9]+\.[0-9]+)','${num2str(str2num($1),''%.3e'')}'),...
        '(e[+-][0-9]+)',' \\times 10^{${strtok($1,''e'')}}'),...
        '((?<=\{[+-])[0-9])','${regexprep($1,''^0*'','''')}')
end