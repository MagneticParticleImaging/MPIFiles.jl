export try_chmod

function try_chmod(path, mode; recursive=true)
  try
    chmod(path,mode,recursive=recursive)
  catch
  end
  return
end

function ishidden(filename::String)
  @static if Sys.isunix()
    s = basename(filename)
    return (!isempty(s) && s[1] == '.')
  else
    # The windows API only supports unicode via UTF-16
    # We therefore have to translate our UTF8 string to UTF16
    # in case it contains non-ascii characters
    utf16 = transcode(UInt16, filename)
    push!(utf16, 0) # Transcode is not null terminated
    # GetFilesAttributesW expects LPCWSTR which is a pointer to 16-bit values
    lpcwstr = pointer(utf16)
    attr = GC.@preserve ccall((:GetFileAttributesW), Cint, (Ptr{Cwchar_t},),lpcwstr)
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
