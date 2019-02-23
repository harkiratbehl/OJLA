function filename = replace_file_ext(filename, ext)
    ind = strfind(filename, '.');
    str = filename(1:ind-1);
    filename = [str '.' ext];

end