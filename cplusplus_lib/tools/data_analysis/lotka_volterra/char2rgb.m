function rgbvec = char2rgb (charcolor)
%
%converts a character color (one of 'r','g','b','c','m','y','k','w') to a 3
%value RGB vector
if (~ischar(charcolor))
    warning('RGB2VEC:NOTC', 'You must pass a character (rgbcmykw)');
    rgbvec = [0 0 0];
    return;
end
rgbvec = zeros(length(charcolor), 3);
for j = 1:length(charcolor)
    switch(lower(charcolor(j)))
        case 'r'
            rgbvec(j,:) = [1 0 0];
        case 'g'
            rgbvec(j,:) = [0 1 0];
        case 'b'
            rgbvec(j,:) = [0 0 1];
        case 'c'
            rgbvec(j,:) = [0 1 1];
        case 'm'
            rgbvec(j,:) = [1 0 1];
        case 'y'
            rgbvec(j,:) = [1 1 0];
        case 'w'
            rgbvec(j,:) = [1 1 1];
        case 'k'
            rgbvec(j,:) = [0 0 0];
        otherwise

    end
end