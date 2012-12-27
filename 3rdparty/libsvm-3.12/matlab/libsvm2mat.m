function params = libsvm2mat(fname)
%LIBSVM2MAT LibSVM to MATLAB conversion
% usage:
%        params = libsvm2mat(model_file_name)
% or if you have a data file you would like to load
%        params = libsvm2mat(data_file_name)

% Copyleft 2003 Thomas Philip Runarsson (26/3)

fid = fopen(fname,'r');
if (fid == -1), error(sprintf('%s not found!',fname')); end
while 1,
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  I = [findstr(tline,' ') length(tline)+1];  
  switch tline(1:I(1)-1),
    case 'svm_type',
      params.svm_type = tline((I(1)+1):(I(2)-1));
    case 'kernel_type', 
      params.kernel_type = tline((I(1)+1):(I(2)-1));
    case 'degree'
      params.degree = sscanf(tline(I(1):end),'%f');
    case 'gamma'
      params.gamma = sscanf(tline(I(1):end),'%f');
    case 'coef0'
      params.coef0 = sscanf(tline(I(1):end),'%f');
    case 'nr_class',
      params.nr_class = sscanf(tline(I(1):end),'%f');
    case 'nr_sv',
      params.nr_sv = sscanf(tline(I(1):end),'%f');
    case 'rho',
      params.rho = sscanf(tline(I(1):end),'%f');
      bias = params.rho;
    case 'nr_sv',
      params.nr_sv = sscanf(tline(I(1):end),'%f');
    case 'total_sv',
      params.total_sv = sscanf(tline(I(1):end),'%f');
    case 'label',
      params.label = sscanf(tline(I(1):end),'%f');
    case 'SV',
      k = 0;
      while 1,
        k = k + 1;
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        tline(find(tline == ':')) = ' ';
        values = sscanf(tline, '%f')';
        params.alpha(k,:) = values(1:params.nr_class-1);
        params.SV(k, values(params.nr_class:2:end)) = values((params.nr_class+1):2:end);
      end
    otherwise, % assume this to be a data file, double check if its true:
      if isempty(str2num(tline(1:I(1)-1))), 
        disp([tline(1:I(1)-1) ' ignored!']);
      else
        k = 0;
        while 1,
          k = k + 1;
          if ~ischar(tline), break, end
          tline(find(tline == ':')) = ' ';
          values = sscanf(tline, '%f')';
          params.y(k,1) = values(1);
          params.X(k, values(2:2:end)) = values(3:2:end);
          tline = fgetl(fid);    
        end    
      end
  end
end
fclose(fid);
