export try_chmod

function try_chmod(path, mode; recursive=true)
  try
    chmod(path,mode,recursive=recursive)
  catch
  end
  return
end

function ishidden(filename::AbstractString)
  @static if Sys.isunix()
    s = basename(filename)
    return (!isempty(s) && s[1] == '.')
  else
    attr = ccall((:GetFileAttributesA), stdcall, Cint, (Ptr{UInt8},),unsafe_string(filename))
    return attr & 0x2 > 0
  end
end


function getNewNumInFolder(path)
  if !isdir(path)
    mkpath(path)
    try_chmod(path, 0o777, recursive=true)
    return 1
  end

  files = readdir(path)
  num = 1
  if length(files) > 0
    for i=1:length(files)
      pref, ext = splitext(files[i])
      num_ = tryparse(Int64, pref)
      if num_ != nothing && num_+1>num
        num = num_+1
      end
    end
  end

  return num
end
