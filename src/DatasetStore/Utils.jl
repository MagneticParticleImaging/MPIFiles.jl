
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
    attr = ccall((:GetFileAttributesA), stdcall, Cint, (Ptr{UInt8},), Base.cconvert(Cstring, filename))
    return attr & 0x2 > 0
  end
end